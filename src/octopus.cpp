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
#include "vcf.hpp"
#include "vcf_utils.hpp"
#include "variant_caller.hpp"

#include "test_common.hpp"
#include "mappable_debug.hpp"

namespace Octopus
{
    size_t count_reads(ReadManager& read_manager, ReferenceGenome& reference)
    {
        size_t result {};
        for (const auto& contig : reference.get_contig_names()) {
            result += read_manager.count_reads(reference.get_contig_region(contig));
        }
        return result;
    }
    
    auto get_contigs(const SearchRegions& regions)
    {
        std::vector<GenomicRegion::ContigNameType> result {};
        result.reserve(regions.size());
        std::transform(std::cbegin(regions), std::cend(regions), std::back_inserter(result),
                       [] (const auto& p) { return p.first; });
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
    
    size_t approx_num_reads(size_t bytes_available)
    {
        return bytes_available / sizeof(AlignedRead);
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
    
    struct GenomeCallingComponents
    {
        ReferenceGenome reference;
        ReadManager read_manager;
        std::vector<SampleIdType> samples;
        SearchRegions regions;
        ReadPipe read_pipe;
        CandidateGeneratorBuilder candidate_generator_builder;
        VcfWriter output;
        
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
        read_pipe {make_read_pipe(this->read_manager, this->samples, options)},
        candidate_generator_builder {Options::make_candidate_generator_builder(options, this->reference)},
        output {std::move(output)}
        {}
        
        ~GenomeCallingComponents() = default;
        
        GenomeCallingComponents(const GenomeCallingComponents&)            = delete;
        GenomeCallingComponents& operator=(const GenomeCallingComponents&) = delete;
        
        GenomeCallingComponents(GenomeCallingComponents&& other) noexcept
        :
        reference {std::move(other.reference)},
        read_manager {std::move(other.read_manager)},
        samples {std::move(other.samples)},
        regions {std::move(other.regions)},
        read_pipe {std::move(other.read_pipe)},
        candidate_generator_builder {std::move(other.candidate_generator_builder)},
        output {std::move(other.output)}
        {
            // need to update or will be pointing to dangling reference
            read_pipe.set_read_manager(read_manager);
            candidate_generator_builder.set_reference(reference);
        }
        
        GenomeCallingComponents& operator=(GenomeCallingComponents&& other) noexcept
        {
            using std::swap;
            swap(reference, other.reference);
            swap(read_manager, other.read_manager);
            swap(samples, other.samples);
            swap(regions, other.regions);
            swap(read_pipe, other.read_pipe);
            swap(candidate_generator_builder, other.candidate_generator_builder);
            swap(output, other.output);
            
            // need to update or will be pointing to dangling reference
            read_pipe.set_read_manager(read_manager);
            candidate_generator_builder.set_reference(reference);
            
            return *this;
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
        
        GenomeCallingComponents result {
            std::move(*reference), std::move(*read_manager),
            std::move(output), options
        };
        
        if (!are_components_valid(result)) return boost::none;
        
        return boost::optional<GenomeCallingComponents> {std::move(result)};
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
        
        ContigCallingComponents(const GenomicRegion::ContigNameType& contig,
                                GenomeCallingComponents& genome_components,
                                const po::variables_map& options)
        :
        reference {genome_components.reference},
        read_manager {genome_components.read_manager},
        regions {genome_components.regions.at(contig)},
        samples {genome_components.samples},
        caller {Options::make_variant_caller(genome_components.reference,
                                             genome_components.read_pipe,
                                             genome_components.candidate_generator_builder,
                                             contig,
                                             options)},
        output {genome_components.output}
        {}
        
        ~ContigCallingComponents() = default;
        
        ContigCallingComponents(const ContigCallingComponents&)            = delete;
        ContigCallingComponents& operator=(const ContigCallingComponents&) = delete;
        ContigCallingComponents(ContigCallingComponents&&)                 = default;
        ContigCallingComponents& operator=(ContigCallingComponents&&)      = default;
    };
    
    void write_final_output_header(GenomeCallingComponents& components)
    {
        components.output.write(make_header(components.samples, get_contigs(components.regions),
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
    
    VcfWriter create_temp_output_file(const GenomicRegion::ContigNameType& contig,
                                      const po::variables_map& options)
    {
        const auto temp_directory = Options::get_temp_file_directory(options);
        
        VcfWriter result {*temp_directory};
        
        return result;
    }
    
    void run_octopus(const po::variables_map& options)
    {
        auto components = collate_genome_calling_components(options);
        
        if (!components) return;
        
        print_startup_info(*components);
        
        write_final_output_header(*components);
        
        for (const auto& p : components->regions) {
            run_octopus_on_contig({p.first, *components, options});
        }
        
        print_final_info(*components);
    }
    
} // namespace Octopus


