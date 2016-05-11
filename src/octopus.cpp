//
//  octopus.cpp
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "octopus.hpp"

#include <vector>
#include <deque>
#include <queue>
#include <map>
#include <iostream>
#include <memory>
#include <algorithm>
#include <numeric>
#include <cstddef>
#include <thread>
#include <future>
#include <functional>
#include <random>
#include <sstream>
#include <chrono>

#include <boost/optional.hpp>

#include "common.hpp"
#include "program_options.hpp"
#include "genomic_region.hpp"
#include "mappable_flat_multi_set.hpp"
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
#include "maths.hpp"
#include "progress_meter.hpp"

#include <cassert>
#include "timers.hpp" // BENCHMARK

#include "logging.hpp"
#include "timing.hpp"

namespace Octopus
{
    void log_startup()
    {
        Logging::InfoLogger log {};
        log << "------------------------------------------------------------------------";
        if (TRACE_MODE) {
            stream(log) << "Octopus v" << Octopus_version << " (trace mode)";
        } else if (DEBUG_MODE) {
            stream(log) << "Octopus v" << Octopus_version << " (debug mode)";
        } else {
            stream(log) << "Octopus v" << Octopus_version;
        }
        log << "Copyright (c) 2016 University of Oxford";
        log << "------------------------------------------------------------------------";
    }
    
    namespace
    {
    template <typename T>
    std::size_t index_of(const std::vector<T>& elements, const T& value)
    {
        return std::distance(std::cbegin(elements), std::find(std::cbegin(elements), std::cend(elements), value));
    }
    
    auto get_contigs(const InputRegionMap& regions, const ReferenceGenome& reference,
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
    
    template <typename Container>
    bool is_in_file_samples(const SampleIdType& sample, const Container& file_samples)
    {
        return std::find(std::cbegin(file_samples), std::cend(file_samples), sample)
                    != std::cend(file_samples);
    }
    
    std::vector<SampleIdType> extract_samples(const po::variables_map& options,
                                              const ReadManager& read_manager)
    {
        auto user_samples = Options::get_user_samples(options);
        auto file_samples = read_manager.get_samples();
        
        if (user_samples) {
            const auto it = std::partition(std::begin(*user_samples), std::end(*user_samples),
                                           [&file_samples] (const auto& sample) {
                                               return is_in_file_samples(sample, file_samples);
                                           });
            
            if (it != std::end(*user_samples)) {
                std::ostringstream ss {};
                ss << "User samples not found in read files: ";
                std::copy(it, std::end(*user_samples), std::ostream_iterator<SampleIdType>(ss, " "));
                Logging::WarningLogger log {};
                log << ss.str();
                user_samples->erase(it, std::end(*user_samples));
            }
            
            return *user_samples;
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
        
        // TODO: find a better place to do this
        vcf_header_builder.add_filter("MODEL", "The caller specific model filter failed");
        vcf_header_builder.add_format("SCR", "2", "Float", "99% credible region of the somatic allele frequency");
        
        return vcf_header_builder.build_once();
    }
    
    VcfHeader make_header(const std::vector<SampleIdType>& samples,
                          const GenomicRegion::ContigNameType& contig,
                          const ReferenceGenome& reference)
    {
        return make_header(samples, std::vector<GenomicRegion::ContigNameType> {contig}, reference);
    }
    
    GenomicRegion::SizeType calculate_total_search_size(const InputRegionMap& regions)
    {
        return std::accumulate(std::cbegin(regions), std::cend(regions), GenomicRegion::SizeType {},
                               [] (const auto curr, const auto& p) {
                                   return curr + sum_region_sizes(p.second);
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
    
    template <typename InputIt, typename RandomGenerator>
    InputIt random_select(InputIt first, InputIt last, RandomGenerator& g) {
        std::uniform_int_distribution<std::size_t> dis(0, std::distance(first, last) - 1);
        std::advance(first, dis(g));
        return first;
    }
    
    template <typename InputIt>
    InputIt random_select(InputIt first, InputIt last) {
        static std::random_device rd {};
        static std::mt19937 gen(rd());
        return random_select(first, last, gen);
    }
    
    auto estimate_read_size(const AlignedRead& read)
    {
        return sizeof(AlignedRead)
            // Now the dynamically allocated bits
            + sequence_size(read) * sizeof(char)
            + sequence_size(read) * sizeof(AlignedRead::QualityType)
            + read.get_cigar_string().size() * sizeof(CigarOperation)
            + contig_name(read).size()
            + (read.is_chimeric() ? sizeof(AlignedRead::NextSegment) : 0);
    }
    
    auto estimate_mean_read_size(const std::vector<SampleIdType>& samples,
                                 const InputRegionMap& input_regions,
                                 ReadManager& read_manager,
                                 unsigned max_sample_size = 1000)
    {
        assert(!input_regions.empty());
        
        const auto num_samples_per_sample = max_sample_size / samples.size();
        
        std::vector<std::size_t> read_sizes {};
        read_sizes.reserve(max_sample_size);
        
        for (const auto& sample : samples) {
            const auto it  = random_select(std::cbegin(input_regions), std::cend(input_regions));
            
            assert(!it->second.empty());
            
            const auto it2 = random_select(std::cbegin(it->second), std::cend(it->second));
            
            const auto test_region = read_manager.find_covered_subregion(sample, *it2, num_samples_per_sample);
            
            const auto reads = read_manager.fetch_reads(sample, test_region);
            
            std::transform(std::cbegin(reads), std::cend(reads), std::back_inserter(read_sizes),
                           estimate_read_size);
        }
        
        if (read_sizes.empty()) {
            return sizeof(AlignedRead) + 300;
        }
        
        return static_cast<std::size_t>(Maths::mean(read_sizes) + Maths::stdev(read_sizes));
    }
    
    std::size_t calculate_max_num_reads(const std::size_t max_buffer_size_in_bytes,
                                        const std::vector<SampleIdType>& samples,
                                        const InputRegionMap& input_regions,
                                        ReadManager& read_manager)
    {
        if (samples.empty()) return 0;
        return max_buffer_size_in_bytes / estimate_mean_read_size(samples, input_regions, read_manager);
    }
    
    class GenomeCallingComponents
    {
    public:
        using Samples = std::vector<SampleIdType>;
        using Contigs = std::vector<GenomicRegion::ContigNameType>;
        
        GenomeCallingComponents() = delete;
        
        explicit GenomeCallingComponents(ReferenceGenome&& reference, ReadManager&& read_manager,
                                         VcfWriter&& output, const po::variables_map& options)
        :
        components_ {std::move(reference), std::move(read_manager), std::move(output), options}
        {}
        
        ~GenomeCallingComponents() = default;
        
        GenomeCallingComponents(const GenomeCallingComponents&)            = delete;
        GenomeCallingComponents& operator=(const GenomeCallingComponents&) = delete;
        
        GenomeCallingComponents(GenomeCallingComponents&& other) noexcept
        :
        components_ {std::move(other.components_)}
        {
            update_dependents();
        }
        
        GenomeCallingComponents& operator=(GenomeCallingComponents&& other) = delete;
        
        const ReferenceGenome& get_reference() const noexcept
        {
            return components_.reference;
        }
        
        ReadManager& get_read_manager() noexcept
        {
            return components_.read_manager;
        }
        
        const ReadManager& get_read_manager() const noexcept
        {
            return components_.read_manager;
        }
        
        ReadPipe& get_read_pipe() noexcept
        {
            return components_.read_pipe;
        }
        
        const ReadPipe& get_read_pipe() const noexcept
        {
            return components_.read_pipe;
        }
        
        const Samples& get_samples() const noexcept
        {
            return components_.samples;
        }
        
        const InputRegionMap& get_search_regions() const noexcept
        {
            return components_.regions;
        }
        
        const Contigs& get_contigs_in_output_order() const noexcept
        {
            return components_.contigs_in_output_order;
        }
        
        VcfWriter& get_output() noexcept
        {
            return components_.output;
        }
        
        const VcfWriter& get_output() const noexcept
        {
            return components_.output;
        }
        
        std::size_t get_read_buffer_size() const noexcept
        {
            return components_.read_buffer_size;
        }
        
        const boost::optional<fs::path>& get_temp_directory() const noexcept
        {
            return components_.temp_directory;
        }
        
        boost::optional<unsigned> get_num_threads() const noexcept
        {
            return components_.num_threads;
        }
        
        const CandidateGeneratorBuilder& get_candidate_generator_builder() const noexcept
        {
            return components_.candidate_generator_builder;
        }
        
        const VariantCallerFactory& get_caller_factory() const noexcept
        {
            return components_.variant_caller_factory;
        }
        
    private:
        struct Components
        {
            Components() = delete;
            
            explicit Components(ReferenceGenome&& reference, ReadManager&& read_manager,
                                VcfWriter&& output, const po::variables_map& options)
            :
            reference {std::move(reference)},
            read_manager {std::move(read_manager)},
            samples {extract_samples(options, this->read_manager)},
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
            num_threads {Options::get_num_threads(options)},
            read_buffer_size {},
            temp_directory {((!num_threads || *num_threads > 1) ? Options::create_temp_file_directory(options) : boost::none)}
            {
                if (!(samples.empty() || regions.empty() || !read_manager.good())) {
                    read_buffer_size = calculate_max_num_reads(Options::get_target_read_buffer_size(options),
                                                               this->samples, this->regions,
                                                               this->read_manager);
                }
            }
            
            ~Components() = default;
            
            Components(const Components&)            = delete;
            Components& operator=(const Components&) = delete;
            Components(Components&&)                 = default;
            Components& operator=(Components&&)      = default;
            
            ReferenceGenome reference;
            ReadManager read_manager;
            Samples samples;
            InputRegionMap regions;
            Contigs contigs_in_output_order;
            ReadPipe read_pipe;
            CandidateGeneratorBuilder candidate_generator_builder;
            VariantCallerFactory variant_caller_factory;
            VcfWriter output;
            boost::optional<unsigned> num_threads;
            std::size_t read_buffer_size;
            boost::optional<fs::path> temp_directory;
        };
        
        Components components_;
        
        void update_dependents() noexcept
        {
            components_.read_pipe.set_read_manager(components_.read_manager);
            components_.candidate_generator_builder.set_reference(components_.reference);
            components_.variant_caller_factory.set_reference(components_.reference);
            components_.variant_caller_factory.set_read_pipe(components_.read_pipe);
            components_.variant_caller_factory.set_candidate_generator_builder(components_.candidate_generator_builder);
        }
    };
    
    bool are_components_valid(const GenomeCallingComponents& components)
    {
        Logging::FatalLogger log {};
        
        if (components.get_samples().empty()) {
            log << "Quiting as no samples were detected";
            return false;
        }
        
        if (components.get_search_regions().empty()) {
            log << "Quiting as there are no input regions";
            return false;
        }
        
        if (components.get_candidate_generator_builder().num_generators() == 0) {
            log << "Quiting as there are no candidate generators";
            return false;
        }
        
        return true;
    }
    
    boost::optional<GenomeCallingComponents>
    collate_genome_calling_components(const po::variables_map& options)
    {
        Logging::FatalLogger log {};
        
        auto reference = Options::make_reference(options);
        
        if (!reference) {
            log << "Quiting as could not make reference genome";
            return boost::none;
        }
        
        auto read_manager = Options::make_read_manager(options);
        
        if (!read_manager) {
            log << "Quiting as there are no read files";
            return boost::none;
        }
        
        auto output = Options::make_output_vcf_writer(options);
        
        if (!output.is_open()) {
            log << "Quiting as could not open output file";
            return boost::none;
        }
        
        try {
            GenomeCallingComponents result {
                std::move(*reference),
                std::move(*read_manager),
                std::move(output),
                options
            };
            
            if (!are_components_valid(result)) {
                return boost::none;
            }
            
            return boost::optional<GenomeCallingComponents> {std::move(result)};
        } catch (const fs::filesystem_error& e) {
            return boost::none; // should already have logged this
        } catch (const std::exception& e) {
            stream(log) << "Could not collate options due to error '" << e.what() << "'";
            return boost::none;
        }
    }
    
    struct ContigCallingComponents
    {
        std::reference_wrapper<const ReferenceGenome> reference;
        std::reference_wrapper<ReadManager> read_manager;
        const InputRegionMap::mapped_type regions;
        std::reference_wrapper<const std::vector<SampleIdType>> samples;
        std::unique_ptr<const VariantCaller> caller;
        std::size_t read_buffer_size;
        std::reference_wrapper<VcfWriter> output;
        
        ContigCallingComponents() = delete;
        
        explicit ContigCallingComponents(const GenomicRegion::ContigNameType& contig,
                                         GenomeCallingComponents& genome_components)
        :
        reference {genome_components.get_reference()},
        read_manager {genome_components.get_read_manager()},
        regions {genome_components.get_search_regions().at(contig)},
        samples {genome_components.get_samples()},
        caller {genome_components.get_caller_factory().make(contig)},
        read_buffer_size {genome_components.get_read_buffer_size()},
        output {genome_components.get_output()}
        {}
        
        explicit ContigCallingComponents(const GenomicRegion::ContigNameType& contig,
                                         VcfWriter& output,
                                         GenomeCallingComponents& genome_components)
        :
        reference {genome_components.get_reference()},
        read_manager {genome_components.get_read_manager()},
        regions {genome_components.get_search_regions().at(contig)},
        samples {genome_components.get_samples()},
        caller {genome_components.get_caller_factory().make(contig)},
        read_buffer_size {genome_components.get_read_buffer_size()},
        output {output}
        {}
        
        ~ContigCallingComponents() = default;
        
        ContigCallingComponents(const ContigCallingComponents&)            = delete;
        ContigCallingComponents& operator=(const ContigCallingComponents&) = delete;
        ContigCallingComponents(ContigCallingComponents&&)                 = default;
        ContigCallingComponents& operator=(ContigCallingComponents&&)      = default;
    };
    
    bool region_has_reads(const GenomicRegion& region, ContigCallingComponents& components)
    {
        return components.read_manager.get().count_reads(components.samples.get(), region) > 0;
    }
    
    void write_final_output_header(GenomeCallingComponents& components, const bool sites_only = false)
    {
        if (sites_only) {
            components.get_output().write(make_header({}, components.get_contigs_in_output_order(),
                                                      components.get_reference()));
        } else {
            components.get_output().write(make_header(components.get_samples(),
                                                      components.get_contigs_in_output_order(),
                                                      components.get_reference()));
        }
    }
    
    void log_startup_info(const GenomeCallingComponents& components)
    {
        std::ostringstream ss {};
        
        if (components.get_samples().size() == 1) {
            ss << "Sample is: ";
        } else {
            ss << "Samples are: ";
        }
        
        std::transform(std::cbegin(components.get_samples()), std::cend(components.get_samples()),
                       std::ostream_iterator<std::string> {ss, " "},
                       [] (const auto& sample) -> std::string {
                           return "\"" + sample + "\"";
                       });
        
        auto str = ss.str();
        str.pop_back(); // the extra whitespace
        
        Logging::InfoLogger log {};
        log << str;
        
        stream(log) << "Writing calls to " << components.get_output().path();
    }
    
    void write_calls(VcfWriter& out, std::deque<VcfRecord>&& calls)
    {
        if (DEBUG_MODE) {
            Logging::DebugLogger log {};
            stream(log) << "Writing " << calls.size() << " calls to VCF";
        }
        for (auto&& call : calls) out.write(std::move(call));
    }
    
    auto find_max_window(const ContigCallingComponents& components,
                        const GenomicRegion& remaining_call_region)
    {
        return components.read_manager.get().find_covered_subregion(components.samples,
                                                                    remaining_call_region,
                                                                    components.read_buffer_size);
    }
    
    auto propose_call_subregion(const ContigCallingComponents& components,
                                const GenomicRegion& remaining_call_region)
    {
        if (is_empty(remaining_call_region)) {
            return remaining_call_region;
        }
        
        const auto max_window = find_max_window(components, remaining_call_region);
        
        if (ends_before(remaining_call_region, max_window)) {
            return remaining_call_region;
        }
        
        return max_window;
    }
    
    auto propose_call_subregion(const ContigCallingComponents& components,
                                const GenomicRegion& current_subregion,
                                const GenomicRegion& input_region)
    {
        assert(contains(input_region, current_subregion));
        return propose_call_subregion(components, right_overhang_region(input_region, current_subregion));
    }
    
    void run_octopus_on_contig(ContigCallingComponents&& components)
    {
        #ifdef BENCHMARK
        init_timers();
        #endif
        
        Logging::InfoLogger log {};
        
        std::deque<VcfRecord> calls;
        
        for (const auto& input_region : components.regions) {
            stream(log) << "Processing input region " << input_region;
            
            ProgressMeter progress_meter {input_region};
            
            auto subregion = propose_call_subregion(components, input_region);
            
            if (is_empty(subregion) && !region_has_reads(input_region, components)) {
                Logging::WarningLogger lg {};
                stream(lg) << "No reads found in input region " << input_region;
            }
            
            while (!is_empty(subregion)) {
                if (DEBUG_MODE) {
                    Logging::DebugLogger lg {};
                    stream(lg) << "Processing subregion " << subregion;
                }
                
                try {
                    calls = components.caller->call(subregion, progress_meter);
                } catch(...) {
                    // TODO: which exceptions can we recover from?
                    throw;
                }
                
                try {
                    write_calls(components.output, std::move(calls));
                } catch(...) {
                    // TODO: which exceptions can we recover from?
                    throw;
                }
                
                subregion = propose_call_subregion(components, subregion, input_region);
            }
        }
        
        #ifdef BENCHMARK
        print_caller_timers();
        #endif
    }
    
    void cleanup(GenomeCallingComponents& components) noexcept
    {
        Logging::InfoLogger log {};
        if (components.get_temp_directory()) {
            try {
                const auto num_files_removed = fs::remove_all(*components.get_temp_directory());
                stream(log) << "Removed " << num_files_removed << " temporary files";
            } catch (const std::exception& e) {
                stream(log) << "Cleanup failed with exception: " << e.what();
            }
        }
    }
    
    void run_octopus_single_threaded(GenomeCallingComponents& components)
    {
        for (const auto& contig : components.get_contigs_in_output_order()) {
            run_octopus_on_contig(ContigCallingComponents {contig, components});
        }
    }
    
    VcfWriter create_unique_temp_output_file(const GenomicRegion& region,
                                             const GenomeCallingComponents& components)
    {
        auto path = *components.get_temp_directory();
        
        const auto& contig = region.get_contig_name();
        const auto begin   = std::to_string(region.get_begin());
        const auto end     = std::to_string(region.get_end());
        
        auto file_name = fs::path {contig + "_" + begin + "-" + end + "_temp.bcf"};
        
        path /= file_name;
        
        return VcfWriter {
            path, make_header(components.get_samples(), contig, components.get_reference())
        };
    }
    
    VcfWriter create_unique_temp_output_file(const GenomicRegion::ContigNameType& contig,
                                             const GenomeCallingComponents& components)
    {
        return create_unique_temp_output_file(components.get_reference().get_contig_region(contig),
                                              components);
    }
    
    auto make_temp_writers(const GenomeCallingComponents& components)
    {
        std::unordered_map<GenomicRegion::ContigNameType, VcfWriter> result {};
        result.reserve(components.get_contigs_in_output_order().size());
        
        for (const auto& contig : components.get_contigs_in_output_order()) {
            result.emplace(contig, create_unique_temp_output_file(contig, components));
        }
        
        return result;
    }
    
    struct Task
    {
        enum class ExecutionPolicy { Threaded, VectorThreaded, None };
        
        Task() = delete;
        
        Task(GenomicRegion region, ExecutionPolicy policy = ExecutionPolicy::None)
        : region {std::move(region)}, policy {policy} {};
        
        GenomicRegion region;
        ExecutionPolicy policy;
    };
    
    struct ContigOrder
    {
        using ContigNameType = GenomicRegion::ContigNameType;
        
        template <typename Container>
        ContigOrder(const Container& contigs)
        : contigs_ {std::cbegin(contigs), std::cend(contigs)} {}
        
        bool operator()(const ContigNameType& lhs, const ContigNameType& rhs) const
        {
            const auto it1 = std::find(std::cbegin(contigs_), std::cend(contigs_), lhs);
            const auto it2 = std::find(std::cbegin(contigs_), std::cend(contigs_), rhs);
            return it1 < it2;
        }
        
    private:
        std::vector<ContigNameType> contigs_;
    };
    
    using TaskQueue = std::queue<Task>;
    using ContigTaskMap = std::map<GenomicRegion::ContigNameType, TaskQueue, ContigOrder>;
    
    TaskQueue divide_work_into_tasks(const ContigCallingComponents& components,
                                     const Task::ExecutionPolicy policy)
    {
        TaskQueue result {};
        
        if (components.regions.empty()) return result;
        
        for (const auto& region : components.regions) {
            auto subregion = propose_call_subregion(components, region);
            
            result.emplace(subregion, policy);
            
            while (!is_empty(subregion)) {
                subregion = propose_call_subregion(components, subregion, region);
                result.emplace(subregion, policy);
            }
        }
        
        return result;
    }
    
    Task::ExecutionPolicy make_execution_policy(const GenomeCallingComponents& components)
    {
        if (components.get_num_threads()) {
            return Task::ExecutionPolicy::None;
        }
        return Task::ExecutionPolicy::Threaded;
    }
    
    ContigTaskMap make_tasks(GenomeCallingComponents& components, const unsigned num_threads)
    {
        const auto policy = make_execution_policy(components);
        
        ContigTaskMap result {ContigOrder {components.get_contigs_in_output_order()}};
        
        for (const auto& contig : components.get_contigs_in_output_order()) {
            ContigCallingComponents contig_components {contig, components};
            contig_components.read_buffer_size /= num_threads;
            result.emplace(contig, divide_work_into_tasks(contig_components, policy));
        }
        
        return result;
    }
    
    unsigned calculate_num_task_threads(const GenomeCallingComponents& components)
    {
        if (components.get_num_threads()) {
            return *components.get_num_threads();
        }
        
        const auto num_hardware_threads = std::thread::hardware_concurrency();
        
        if (num_hardware_threads > 0) return num_hardware_threads;
        
        const auto num_files = components.get_read_manager().num_files();
        
        return std::min(num_files, 8u);
    }
    
    Task pop(ContigTaskMap& tasks)
    {
        assert(!tasks.empty());
        auto it = std::begin(tasks);
        auto result = it->second.front();
        it->second.pop();
        if (it->second.empty()) {
            tasks.erase(it->first);
        }
        return result;
    }
    
    template<typename R>
    bool is_ready(const std::future<R>& f)
    {
        assert(f.valid());
        return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
    }
    
    struct CompletedTask
    {
        CompletedTask(Task task, std::deque<VcfRecord>&& calls)
        : task {std::move(task)}, calls {std::move(calls)} {}
        Task task;
        std::deque<VcfRecord> calls;
    };
    
    auto run(const Task& task, GenomeCallingComponents& components, ProgressMeter& progress_meter)
    {
        if (DEBUG_MODE) {
            Logging::DebugLogger log {};
            stream(log) << "Spawning task with region " << task.region;
        }
        
        ContigCallingComponents contig_components {task.region.get_contig_name(), components};
        
        return std::async(std::launch::async,
                          [&task, &progress_meter, components = std::move(contig_components)]
                            () -> CompletedTask {
                              return CompletedTask {
                                  task,
                                  components.caller->call(task.region, progress_meter)
                              };
                          });
    }
    
    template <typename K>
    auto extract(std::unordered_map<K, VcfWriter>& writers)
    {
        std::vector<VcfWriter> result {};
        result.reserve(writers.size());
        
        for (auto& p : writers) result.emplace_back(std::move(p.second));
        
        writers.clear();
        
        return result;
    }
    
    template <typename K>
    auto convert_temp_writers(std::unordered_map<K, VcfWriter>& writers)
    {
        auto tmp = extract(writers);
        return writers_to_readers(tmp);
    }
    
    void run_octopus_multi_threaded(GenomeCallingComponents& components)
    {
        const auto num_task_threads = calculate_num_task_threads(components);
        
        auto tasks = make_tasks(components, num_task_threads);
        
        auto temp_writers = make_temp_writers(components);
        
        std::vector<std::future<CompletedTask>> futures(num_task_threads);
        
        ProgressMeter meter {components.get_search_regions()};
        
        while (!tasks.empty()) {
            for (auto& fut : futures) {
                if (fut.valid()) {
                    if (is_ready(fut)) {
                        try {
                            auto result = fut.get();
                            write_calls(temp_writers.at(result.task.region.get_contig_name()),
                                        std::move(result.calls));
                        } catch (...) {
                            // TODO: which exceptions can we recover from?
                            throw;
                        }
                    }
                } else if (!tasks.empty()) {
                    fut = run(pop(tasks), components, meter);
                }
            }
        }
        
        for (auto& fut : futures) {
            if (fut.valid()) {
                try {
                    auto result = fut.get();
                    write_calls(temp_writers.at(result.task.region.get_contig_name()),
                                std::move(result.calls));
                } catch (...) {
                    // TODO: which exceptions can we recover from?
                    throw;
                }
            }
        }
        
        auto results = convert_temp_writers(temp_writers);
        
        index_vcfs(results);
        
        merge(results, components.get_contigs_in_output_order(), components.get_output());
    }
    } // namespace
    
    bool is_multithreaded(const GenomeCallingComponents& components)
    {
        return !components.get_num_threads() || *components.get_num_threads() > 1;
    }
    
    void log_final_info(const GenomeCallingComponents& components, const TimeInterval& runtime)
    {
        Logging::InfoLogger log {};
        stream(log) << "Finished processing " << calculate_total_search_size(components.get_search_regions())
                    << "bp, total runtime " << runtime;
    }
    
    void run_octopus(po::variables_map& options)
    {
        DEBUG_MODE = Options::is_debug_mode(options);
        TRACE_MODE = Options::is_trace_mode(options);
        
        log_startup();
        
        const auto start = std::chrono::system_clock::now();
        
        auto components = collate_genome_calling_components(options);
        
        if (!components) return;
        
        auto end = std::chrono::system_clock::now();
        
        Logging::InfoLogger log {};
        
        stream(log) << "Done initialising calling components in " << TimeInterval {start, end};
        
        log_startup_info(*components);
        
        write_final_output_header(*components, Options::call_sites_only(options));
        
        options.clear();
        
        if (is_multithreaded(*components)) {
            run_octopus_multi_threaded(*components);
        } else {
            run_octopus_single_threaded(*components);
        }
        
//        try {
//            log_startup_info(*components);
//            
//            write_final_output_header(*components, Options::call_sites_only(options));
//            
//            options.clear();
//            
//            if (is_multithreaded(*components)) {
//                run_octopus_multi_threaded(*components);
//            } else {
//                run_octopus_single_threaded(*components);
//            }
//        } catch (const std::exception& e) {
//            Logging::FatalLogger lg {};
//            stream(lg) << "Encountered exception '" << e.what() << "'. Attempting to cleanup...";
//            cleanup(*components);
//            if (DEBUG_MODE) {
//                stream(lg) << "Cleanup successful. Please send log file to dcooke@well.ox.ac.uk";
//            } else {
//                stream(lg) << "Cleanup successful. Please re-run in debug mode (option --debug) and send"
//                                " log file to " << Octopus_bug_email;
//            }
//            return;
//        }
        
        cleanup(*components);
        
        end = std::chrono::system_clock::now();
        
        log_final_info(*components, TimeInterval {start, end});
    }
    
} // namespace Octopus


