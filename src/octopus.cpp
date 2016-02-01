//
//  octopus.cpp
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "octopus.hpp"

#include <vector>
#include <iostream>
#include <memory>
#include <algorithm>
#include <numeric>
#include <cstddef>
#include <thread>
#include <future>
#include <functional>

#include <boost/optional.hpp>

#include "common.hpp"
#include "program_options.hpp"
#include "genomic_region.hpp"
#include "mappable_set.hpp"
#include "mappable_map.hpp"
#include "reference_genome.hpp"
#include "read_manager.hpp"
#include "read_filter.hpp"
#include "downsampler.hpp"
#include "read_transform.hpp"
#include "read_pipe.hpp"
#include "read_utils.hpp"
#include "candidate_generator_builder.hpp"
#include "variant_caller_factory.hpp"
#include "variant_caller.hpp"
#include "vcf.hpp"

#include "test_common.hpp"
#include "mappable_debug.hpp"

namespace Octopus
{
    template <typename T>
    static std::size_t index_of(const std::vector<T>& elements, const T& value)
    {
        return std::distance(std::cbegin(elements), std::find(std::cbegin(elements), std::cend(elements), value));
    }
    
    auto get_contigs(const SearchRegions& regions, const ReferenceGenome& reference,
                     const Options::ContigOutputOrder order)
    {
        using Options::ContigOutputOrder;
        
        using ContigNameType = GenomicRegion::ContigNameType;
        
        std::vector<ContigNameType> result {};
        result.reserve(regions.size());
        
        std::transform(std::cbegin(regions), std::cend(regions), std::back_inserter(result),
                       [] (const auto& p) { return p.first; });
        
        switch (order) {
            case ContigOutputOrder::LexicographicalAscending:
                std::sort(std::begin(result), std::end(result));
                break;
            case ContigOutputOrder::LexicographicalDescending:
                std::sort(std::begin(result), std::end(result), std::greater<ContigNameType>());
                break;
            case ContigOutputOrder::ContigSizeAscending:
                std::sort(std::begin(result), std::end(result),
                          [&] (const auto& lhs, const auto& rhs) {
                              return reference.get_contig_size(lhs) < reference.get_contig_size(rhs);
                          });
                break;
            case ContigOutputOrder::ContigSizeDescending:
                std::sort(std::begin(result), std::end(result),
                          [&] (const auto& lhs, const auto& rhs) {
                              return reference.get_contig_size(lhs) > reference.get_contig_size(rhs);
                          });
                break;
            case ContigOutputOrder::AsInReferenceIndex:
            {
                const auto reference_contigs = reference.get_contig_names();
                std::sort(std::begin(result), std::end(result),
                          [&] (const auto& lhs, const auto& rhs) {
                              return index_of(reference_contigs, lhs) < index_of(reference_contigs, rhs);
                          });
                break;
            }
            case ContigOutputOrder::AsInReferenceIndexReversed:
            {
                const auto reference_contigs = reference.get_contig_names();
                std::sort(std::begin(result), std::end(result),
                          [&] (const auto& lhs, const auto& rhs) {
                              return index_of(reference_contigs, lhs) > index_of(reference_contigs, rhs);
                          });
                break;
            }
            case ContigOutputOrder::Unspecified:
                break;
        }
        
        return result;
    }
    
    std::vector<SampleIdType> get_samples(const po::variables_map& options,
                                          const ReadManager& read_manager)
    {
        auto user_samples = Options::get_samples(options);
        auto file_samples = read_manager.get_samples();
        
        if (!user_samples.empty()) {
            std::vector<SampleIdType> bad_samples {};
            
            for (const auto& user_sample : user_samples) {
                if (std::find(std::cbegin(file_samples), std::cend(file_samples),
                              user_sample) == std::cend(file_samples)) {
                    bad_samples.push_back(user_sample);
                }
            }
            
            if (!bad_samples.empty()) {
                std::cout << "Octopus: input samples not present in read files: ";
                std::copy(std::cbegin(bad_samples), std::cend(bad_samples),
                          std::ostream_iterator<SampleIdType>(std::cout, " "));
                std::cout << std::endl;
                return {};
            }
            
            return user_samples;
        } else {
            return file_samples;
        }
    }
    
    VcfHeader make_header(const std::vector<SampleIdType>& samples,
                          const std::vector<GenomicRegion::ContigNameType>& contigs,
                          const ReferenceGenome& reference)
    {
        auto vcf_header_builder = get_default_header_builder().set_samples(samples);
        for (const auto& contig : contigs) {
            vcf_header_builder.add_contig(contig, {{"length", std::to_string(reference.get_contig_size(contig))}});
        }
        vcf_header_builder.add_basic_field("reference", reference.get_name());
        vcf_header_builder.add_structured_field("Octopus", {{"some", "option"}});
        
        return vcf_header_builder.build_once();
    }
    
    VcfHeader make_header(const std::vector<SampleIdType>& samples,
                          const GenomicRegion::ContigNameType& contig,
                          const ReferenceGenome& reference)
    {
        return make_header(samples, std::vector<GenomicRegion::ContigNameType> {contig}, reference);
    }
    
    GenomicRegion::SizeType calculate_total_search_size(const SearchRegions& regions)
    {
        return std::accumulate(std::cbegin(regions), std::cend(regions), GenomicRegion::SizeType {},
                               [] (const auto curr, const auto& p) {
                                   return curr + sum_sizes(p.second);
                               });
    }
    
    ReadPipe make_read_pipe(ReadManager& read_manager, std::vector<SampleIdType> samples,
                            const po::variables_map& options)
    {
        return ReadPipe {
            read_manager,
            Options::make_read_filter(options),
            Options::make_downsampler(options),
            Options::make_read_transform(options),
            std::move(samples)
        };
    }
    
    class GenomeCallingComponents
    {
    public:
        ReferenceGenome reference;
        ReadManager read_manager;
        std::vector<SampleIdType> samples;
        SearchRegions regions;
        std::vector<GenomicRegion::ContigNameType> contigs_in_output_order;
        ReadPipe read_pipe;
        CandidateGeneratorBuilder candidate_generator_builder;
        VariantCallerFactory variant_caller_factory;
        VcfWriter output;
        bool is_threaded;
        boost::optional<fs::path> temp_directory;
        
        GenomeCallingComponents() = delete;
        
        explicit GenomeCallingComponents(ReferenceGenome&& reference,
                                         ReadManager&& read_manager,
                                         VcfWriter&& output,
                                         const po::variables_map& options)
        :
        reference {std::move(reference)},
        read_manager {std::move(read_manager)},
        samples {get_samples(options, this->read_manager)},
        regions {Options::get_search_regions(options, this->reference)},
        contigs_in_output_order {get_contigs(this->regions, this->reference,
                                             Options::get_contig_output_order(options))},
        read_pipe {make_read_pipe(this->read_manager, this->samples, options)},
        candidate_generator_builder {Options::make_candidate_generator_builder(options, this->reference)},
        variant_caller_factory {Options::make_variant_caller_factory(this->reference,
                                                                     this->read_pipe,
                                                                     this->candidate_generator_builder,
                                                                     this->regions,
                                                                     options)},
        output {std::move(output)},
        is_threaded {Options::is_threading_allowed(options)},
        temp_directory {is_threaded ? Options::create_temp_file_directory(options) : boost::none}
        {}
        
        ~GenomeCallingComponents() = default;
        
        GenomeCallingComponents(const GenomeCallingComponents&)            = delete;
        GenomeCallingComponents& operator=(const GenomeCallingComponents&) = delete;
        
        GenomeCallingComponents(GenomeCallingComponents&& other) noexcept
        :
        reference                   {std::move(other.reference)},
        read_manager                {std::move(other.read_manager)},
        samples                     {std::move(other.samples)},
        regions                     {std::move(other.regions)},
        contigs_in_output_order     {std::move(other.contigs_in_output_order)},
        read_pipe                   {std::move(other.read_pipe)},
        candidate_generator_builder {std::move(other.candidate_generator_builder)},
        variant_caller_factory      {std::move(other.variant_caller_factory)},
        output                      {std::move(other.output)},
        is_threaded                 {std::move(other.is_threaded)},
        temp_directory              {std::move(other.temp_directory)}
        {
            update_dependents();
        }
        
        GenomeCallingComponents& operator=(GenomeCallingComponents&& other) noexcept
        {
            using std::swap;
            swap(reference,                   other.reference);
            swap(read_manager,                other.read_manager);
            swap(samples,                     other.samples);
            swap(regions,                     other.regions);
            swap(contigs_in_output_order,     other.contigs_in_output_order);
            swap(read_pipe,                   other.read_pipe);
            swap(candidate_generator_builder, other.candidate_generator_builder);
            swap(variant_caller_factory,      other.variant_caller_factory);
            swap(output,                      other.output);
            swap(is_threaded,                 other.is_threaded);
            swap(temp_directory,              other.temp_directory);
            update_dependents();
            return *this;
        }
        
    private:
        void update_dependents() noexcept
        {
            read_pipe.set_read_manager(read_manager);
            candidate_generator_builder.set_reference(reference);
            variant_caller_factory.set_reference(reference);
            variant_caller_factory.set_read_pipe(read_pipe);
            variant_caller_factory.set_candidate_generator_builder(candidate_generator_builder);
        }
    };
    
    bool are_components_valid(const GenomeCallingComponents& components)
    {
        using std::cout; using std::endl;
        
        if (components.samples.empty()) {
            cout << "Octopus: quiting as no samples were found" << endl;
            return false;
        }
        
        if (components.regions.empty()) {
            cout << "Octopus: quiting as got no input regions" << endl;
            return false;
        }
        
        if (components.candidate_generator_builder.num_generators() == 0) {
            cout << "Octopus: quiting as there are no candidate generators" << endl;
            return false;
        }
        
        return true;
    }
    
    boost::optional<GenomeCallingComponents> collate_genome_calling_components(const po::variables_map& options)
    {
        using std::cout; using std::endl;
        
        auto reference = Options::make_reference(options);
        
        if (!reference) {
            cout << "Octopus: quiting as could not make reference genome" << endl;
            return boost::none;
        }
        
        auto read_manager = Options::make_read_manager(options);
        
        if (!read_manager) {
            cout << "Octopus: quiting as could not load read files" << endl;
            return boost::none;
        }
        
        auto output = Options::make_output_vcf_writer(options);
        
        if (!output.is_open()) {
            cout << "Octopus: quiting as could not open output file" << endl;
            return boost::none;
        }
        
        try {
            GenomeCallingComponents result {
                std::move(*reference),
                std::move(*read_manager),
                std::move(output),
                options
            };
            
            if (!are_components_valid(result)) return boost::none;
            
            return boost::optional<GenomeCallingComponents> {std::move(result)};
        } catch (...) {
            std::cout << "Octopus: could not collate options" << std::endl;
            return boost::none;
        }
    }
    
    struct ContigCallingComponents
    {
        std::reference_wrapper<const ReferenceGenome> reference;
        std::reference_wrapper<ReadManager> read_manager;
        const MappableSet<GenomicRegion> regions;
        std::reference_wrapper<const std::vector<SampleIdType>> samples;
        std::unique_ptr<const VariantCaller> caller;
        std::reference_wrapper<VcfWriter> output;
        
        ContigCallingComponents() = delete;
        
        explicit ContigCallingComponents(const GenomicRegion::ContigNameType& contig,
                                         GenomeCallingComponents& genome_components)
        :
        reference {genome_components.reference},
        read_manager {genome_components.read_manager},
        regions {genome_components.regions.at(contig)},
        samples {genome_components.samples},
        caller {genome_components.variant_caller_factory.make(contig)},
        output {genome_components.output}
        {}
        
        explicit ContigCallingComponents(const GenomicRegion::ContigNameType& contig,
                                         VcfWriter& output,
                                         GenomeCallingComponents& genome_components)
        :
        reference {genome_components.reference},
        read_manager {genome_components.read_manager},
        regions {genome_components.regions.at(contig)},
        samples {genome_components.samples},
        caller {genome_components.variant_caller_factory.make(contig)},
        output {output}
        {}
        
        ~ContigCallingComponents() = default;
        
        ContigCallingComponents(const ContigCallingComponents&)            = delete;
        ContigCallingComponents& operator=(const ContigCallingComponents&) = delete;
        ContigCallingComponents(ContigCallingComponents&&)                 = default;
        ContigCallingComponents& operator=(ContigCallingComponents&&)      = default;
    };
    
    void write_final_output_header(GenomeCallingComponents& components)
    {
        components.output.write(make_header(components.samples, components.contigs_in_output_order,
                                            components.reference));
    }
    
    void print_startup_info(const GenomeCallingComponents& components)
    {
        using std::cout; using std::endl;
        cout << "Octopus: calling variants in " << components.samples.size() << " samples" << endl;
        cout << "Octopus: writing calls to " << components.output.path() << endl;
    }
    
    void write_calls(VcfWriter& out, std::vector<VcfRecord>&& calls)
    {
        std::cout << "Octopus: writing " << calls.size() << " calls to VCF" << std::endl;
        for (auto&& call : calls) out.write(std::move(call));
    }
    
    void run_octopus_on_contig(ContigCallingComponents&& components)
    {
        using std::cout; using std::endl;
        
        const size_t max_reads = 1'000'000;
        
        size_t num_buffered_reads {};
        
        for (const auto& region : components.regions) {
            cout << "Octopus: processing input region " << region << endl;
            
            auto subregion = components.read_manager.get().find_covered_subregion(components.samples,
                                                                                  region, max_reads);
            
            while (get_begin(subregion) != get_end(region)) {
                cout << "Octopus: processing subregion " << subregion << endl;
                
                write_calls(components.output, components.caller->call_variants(subregion));
                
                num_buffered_reads = components.caller->num_buffered_reads();
                
                subregion = get_right_overhang(region, subregion);
                subregion = components.read_manager.get().find_covered_subregion(components.samples,
                                                                                 subregion, max_reads);
            }
        }
    }
    
    void print_final_info(const GenomeCallingComponents& components)
    {
        std::cout << "Octopus: processed " << calculate_total_search_size(components.regions) << "bp" << std::endl;
    }
    
    void cleanup(GenomeCallingComponents& components) noexcept
    {
        using std::cout; using std::endl;
        if (components.temp_directory) {
            try {
                const auto num_files_removed = fs::remove_all(*components.temp_directory);
                cout << "Octopus: removed " << num_files_removed << " temporary files" << endl;
            } catch (std::exception& e) {
                cout << "Octopus: cleanup failed with error " << e.what() << endl;
            }
        }
    }
    
    void run_octopus_single_threaded(GenomeCallingComponents& components)
    {
        for (const auto& contig : components.contigs_in_output_order) {
            run_octopus_on_contig(ContigCallingComponents {contig, components});
        }
    }
    
    VcfWriter create_unique_temp_output_file(const GenomicRegion& region,
                                             const GenomeCallingComponents& components)
    {
        auto path = *components.temp_directory;
        const auto& contig = region.get_contig_name();
        const auto begin   = std::to_string(region.get_begin());
        const auto end     = std::to_string(region.get_end());
        
        auto file_name = fs::path {contig + "_" + begin + "-" + end + "_temp.bcf"};
        
        path /= file_name;
        
        return VcfWriter {path, make_header(components.samples, contig, components.reference)};
    }
    
    VcfWriter create_unique_temp_output_file(const GenomicRegion::ContigNameType& contig,
                                             const GenomeCallingComponents& components)
    {
        return create_unique_temp_output_file(components.reference.get_contig_region(contig), components);
    }
    
    std::vector<VcfWriter> create_temp_writers(const GenomeCallingComponents& components)
    {
        std::vector<VcfWriter> result {};
        result.reserve(components.contigs_in_output_order.size());
        for (const auto& contig : components.contigs_in_output_order) {
            result.emplace_back(create_unique_temp_output_file(contig, components));
        }
        return result;
    }
    
    void run_octopus_multi_threaded(GenomeCallingComponents& components)
    {
        auto temp_writers = create_temp_writers(components);
        
        std::vector<std::future<void>> tasks {};
        
        const auto& contigs = components.contigs_in_output_order;
        
        std::transform(std::cbegin(contigs), std::cend(contigs), std::begin(temp_writers),
                       std::back_inserter(tasks),
                       [&] (const auto& contig, auto& writer) {
                           return std::async(run_octopus_on_contig,
                                             ContigCallingComponents {contig, writer, components});
                       });
        
        for (auto& task : tasks) task.get();
        
        auto results = writers_to_readers(temp_writers);
        
        index_vcfs(results);
        
        merge(results, contigs, components.output);
    }
    
    void run_octopus(po::variables_map& options)
    {
        auto components = collate_genome_calling_components(options);
        
        if (!components) return;
        
        options.clear();
        
        print_startup_info(*components);
        
        write_final_output_header(*components);
        
        if (components->is_threaded) {
            run_octopus_multi_threaded(*components);
        } else {
            run_octopus_single_threaded(*components);
        }
        
//        try {
//            print_startup_info(*components);
//            
//            write_final_output_header(*components);
//            
//            if (components->is_threaded) {
//                run_octopus_multi_threaded(*components);
//            } else {
//                run_octopus_single_threaded(*components);
//            }
//        } catch (...) {
//            std::cout << "Octopus: encountered fatal error. Attempting to cleanup..." << std::endl;
//            cleanup(*components);
//            throw;
//        }
        
        print_final_info(*components);
        cleanup(*components);
    }
    
} // namespace Octopus


