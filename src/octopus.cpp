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
#include <condition_variable>
#include <mutex>
#include <atomic>
#include <array>
#include <set>
#include <typeinfo>
#include <cassert>

#include <boost/optional.hpp>

#include "common.hpp"
#include "program_options.hpp"
#include "genomic_region.hpp"
#include "mappable.hpp"
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
#include "vcf_header_factory.hpp"
#include "octopus_vcf.hpp"
#include "variant_caller_factory.hpp"
#include "variant_caller.hpp"
#include "vcf.hpp"
#include "maths.hpp"
#include "progress_meter.hpp"
#include "logging.hpp"
#include "timing.hpp"

#include "variant_call_filter.hpp"
#include "read_transformations.hpp"

#include "timers.hpp" // BENCHMARK

#include "genotype_reader.hpp"

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
                              return reference.contig_size(lhs) < reference.contig_size(rhs);
                          });
                break;
            case ContigOutputOrder::ContigSizeDescending:
                std::sort(std::begin(result), std::end(result),
                          [&] (const auto& lhs, const auto& rhs) {
                              return reference.contig_size(lhs) > reference.contig_size(rhs);
                          });
                break;
            case ContigOutputOrder::AsInReferenceIndex:
            {
                const auto reference_contigs = reference.contig_names();
                std::sort(std::begin(result), std::end(result),
                          [&] (const auto& lhs, const auto& rhs) {
                              return index_of(reference_contigs, lhs) < index_of(reference_contigs, rhs);
                          });
                break;
            }
            case ContigOutputOrder::AsInReferenceIndexReversed:
            {
                const auto reference_contigs = reference.contig_names();
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
        auto file_samples = read_manager.samples();
        
        if (user_samples) {
            const auto it = std::partition(std::begin(*user_samples), std::end(*user_samples),
                                           [&file_samples] (const auto& sample) {
                                               return is_in_file_samples(sample, file_samples);
                                           });
            
            const auto num_not_found = std::distance(it, std::end(*user_samples));
            
            if (num_not_found > 0) {
                std::ostringstream ss {};
                ss << "The requested calling sample";
                if (num_not_found > 1) ss << 's';
                ss << " ";
                std::transform(it, std::end(*user_samples), std::ostream_iterator<SampleIdType>(ss, ", "),
                               [] (auto sample) { return "'" + sample + "'"; });
                if (num_not_found == 1) {
                    ss << "is";
                } else {
                    ss << "are";
                }
                ss << " not present in any of the read files";
                Logging::WarningLogger log {};
                log << ss.str();
                user_samples->erase(it, std::end(*user_samples));
            }
            
            return *user_samples;
        } else {
            return file_samples;
        }
    }
    
    using CallTypeSet = std::set<std::type_index>;
    
    VcfHeader make_vcf_header(const std::vector<SampleIdType>& samples,
                              const std::vector<GenomicRegion::ContigNameType>& contigs,
                              const ReferenceGenome& reference,
                              const CallTypeSet& call_types)
    {
        auto builder = Vcf::make_octopus_header_template().set_samples(samples);
        
        for (const auto& contig : contigs) {
            builder.add_contig(contig, {{"length", std::to_string(reference.contig_size(contig))}});
        }
        
        builder.add_basic_field("reference", reference.name());
        builder.add_structured_field("Octopus", {{"some", "option"}});
        
        VcfHeaderFactory factory {};
        
        for (const auto& type : call_types) {
            factory.register_call_type(type);
        }
        
        factory.annotate(builder);
        
        return builder.build_once();
    }
    
    VcfHeader make_vcf_header(const std::vector<SampleIdType>& samples,
                              const GenomicRegion::ContigNameType& contig,
                              const ReferenceGenome& reference,
                              const CallTypeSet& call_types)
    {
        return make_vcf_header(samples, std::vector<GenomicRegion::ContigNameType> {contig},
                               reference, call_types);
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
            read_manager, Options::make_read_transform(options),
            Options::make_read_filter(options), Options::make_downsampler(options),
            std::move(samples)
        };
    }
    
    template <typename ForwardIt, typename RandomGenerator>
    ForwardIt random_select(ForwardIt first, ForwardIt last, RandomGenerator& g) {
        std::uniform_int_distribution<std::size_t> dis(0, std::distance(first, last) - 1);
        std::advance(first, dis(g));
        return first;
    }
    
    template <typename ForwardIt>
    ForwardIt random_select(ForwardIt first, ForwardIt last) {
        static std::default_random_engine gen {};
        return random_select(first, last, gen);
    }
    
    auto estimate_read_size(const AlignedRead& read)
    {
        return sizeof(AlignedRead)
            // Now the dynamically allocated bits
            + sequence_size(read) * sizeof(char)
            + sequence_size(read) * sizeof(AlignedRead::QualityType)
            + read.cigar_string().size() * sizeof(CigarOperation)
            + contig_name(read).size()
            + (read.has_other_segment() ? sizeof(AlignedRead::NextSegment) : 0);
    }
    
    auto estimate_mean_read_size(const std::vector<SampleIdType>& samples,
                                 const InputRegionMap& input_regions,
                                 ReadManager& read_manager,
                                 unsigned max_sample_size = 1000)
    {
        assert(!input_regions.empty());
        
        const auto num_samples_per_sample = max_sample_size / samples.size();
        
        std::deque<std::size_t> read_sizes {};
        
        for (const auto& sample : samples) {
            const auto it  = random_select(std::cbegin(input_regions), std::cend(input_regions));
            
            assert(!it->second.empty());
            
            const auto it2 = random_select(std::cbegin(it->second), std::cend(it->second));
            
            auto test_region = read_manager.find_covered_subregion(sample, *it2, num_samples_per_sample);
            
            if (is_empty(test_region)) {
                test_region = expand_rhs(test_region, 1);
            }
            
            const auto reads = read_manager.fetch_reads(sample, test_region);
            
            std::transform(std::cbegin(reads), std::cend(reads), std::back_inserter(read_sizes),
                           estimate_read_size);
        }
        
        if (read_sizes.empty()) {
            Logging::WarningLogger log {};
            log << "Could not estimate read size from data, resorting to default";
            
            return sizeof(AlignedRead) + 300;
        }
        
        return static_cast<std::size_t>(Maths::mean(read_sizes) + Maths::stdev(read_sizes));
    }
    
    std::size_t calculate_max_num_reads(const std::size_t max_buffer_bytes,
                                        const std::vector<SampleIdType>& samples,
                                        const InputRegionMap& input_regions,
                                        ReadManager& read_manager)
    {
        if (samples.empty()) return 0;
        static constexpr std::size_t Min_buffer_bytes {1'000'000};
        return std::max(max_buffer_bytes, Min_buffer_bytes)
            / estimate_mean_read_size(samples, input_regions, read_manager);
    }
    
    class GenomeCallingComponents
    {
    public:
        using Samples = std::vector<SampleIdType>;
        using Contigs = std::vector<GenomicRegion::ContigNameType>;
        
        GenomeCallingComponents() = delete;
        
        GenomeCallingComponents(ReferenceGenome&& reference, ReadManager&& read_manager,
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
        
        const ReferenceGenome& reference() const noexcept
        {
            return components_.reference;
        }
        
        ReadManager& read_manager() noexcept
        {
            return components_.read_manager;
        }
        
        const ReadManager& read_manager() const noexcept
        {
            return components_.read_manager;
        }
        
        ReadPipe& read_pipe() noexcept
        {
            return components_.read_pipe;
        }
        
        const ReadPipe& read_pipe() const noexcept
        {
            return components_.read_pipe;
        }
        
        const Samples& samples() const noexcept
        {
            return components_.samples;
        }
        
        const InputRegionMap& search_regions() const noexcept
        {
            return components_.regions;
        }
        
        const Contigs& contigs_in_output_order() const noexcept
        {
            return components_.contigs_in_output_order;
        }
        
        VcfWriter& output() noexcept
        {
            return components_.output;
        }
        
        const VcfWriter& output() const noexcept
        {
            return components_.output;
        }
        
        std::size_t read_buffer_size() const noexcept
        {
            return components_.read_buffer_size;
        }
        
        const boost::optional<fs::path>& temp_directory() const noexcept
        {
            return components_.temp_directory;
        }
        
        boost::optional<unsigned> num_threads() const noexcept
        {
            return components_.num_threads;
        }
        
        const CandidateGeneratorBuilder& candidate_generator_builder() const noexcept
        {
            return components_.candidate_generator_builder;
        }
        
        const VariantCallerFactory& caller_factory() const noexcept
        {
            return components_.variant_caller_factory;
        }
        
        ProgressMeter& progress_meter() noexcept
        {
            return components_.progress_meter;
        }
        
    private:
        struct Components
        {
            Components() = delete;
            
            Components(ReferenceGenome&& reference, ReadManager&& read_manager,
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
            temp_directory {((!num_threads || *num_threads > 1) ? Options::create_temp_file_directory(options) : boost::none)},
            progress_meter {regions}
            {
                const auto num_bp_to_process = sum_region_sizes(regions);
                
                if (num_bp_to_process < 100000000) {
                    progress_meter.set_percent_block_size(1.0);
                } else if (num_bp_to_process < 1000000000) {
                    progress_meter.set_percent_block_size(0.5);
                } else {
                    progress_meter.set_percent_block_size(0.1);
                }
                
                if (!samples.empty() && !regions.empty() && read_manager.good()) {
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
            ProgressMeter progress_meter;
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
        
        if (components.samples().empty()) {
            log << "No samples detected - at least one is required for calling";
            return false;
        }
        
        if (components.search_regions().empty()) {
            log << "There are no input regions - at least one is required for calling";
            return false;
        }
        
        if (components.candidate_generator_builder().num_generators() == 0) {
            log << "There are no candidate generators - at least one is required for calling";
            return false;
        }
        
        return true;
    }
    
    void cleanup(GenomeCallingComponents& components) noexcept
    {
        Logging::InfoLogger log {};
        if (components.temp_directory()) {
            try {
                const auto num_files_removed = fs::remove_all(*components.temp_directory());
                stream(log) << "Removed " << num_files_removed << " temporary files";
                } catch (const std::exception& e) {
                    stream(log) << "Cleanup failed with exception: " << e.what();
                }
        }
    }
    
    boost::optional<GenomeCallingComponents>
    collate_genome_calling_components(const po::variables_map& options)
    {
        try {
            auto reference = Options::make_reference(options);
            
            auto read_manager = Options::make_read_manager(options);
            
            auto output = Options::make_output_vcf_writer(options);
            
            if (!output.is_open()) {
                return boost::none;
            }
            
            GenomeCallingComponents result {
                std::move(reference),
                std::move(read_manager),
                std::move(output),
                options
            };
            
            if (!are_components_valid(result)) {
                cleanup(result);
                return boost::none;
            }
            
            return boost::optional<GenomeCallingComponents> {std::move(result)};
        } catch (const fs::filesystem_error& e) {
            return boost::none; // should already have logged this
        } catch (const std::exception& e) {
            Logging::FatalLogger log {};
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
        std::reference_wrapper<ProgressMeter> progress_meter;
        
        ContigCallingComponents() = delete;
        
        ContigCallingComponents(const GenomicRegion::ContigNameType& contig,
                                GenomeCallingComponents& genome_components)
        :
        reference {genome_components.reference()},
        read_manager {genome_components.read_manager()},
        regions {genome_components.search_regions().at(contig)},
        samples {genome_components.samples()},
        caller {genome_components.caller_factory().make(contig)},
        read_buffer_size {genome_components.read_buffer_size()},
        output {genome_components.output()},
        progress_meter {genome_components.progress_meter()}
        {}
        
        ContigCallingComponents(const GenomicRegion::ContigNameType& contig, VcfWriter& output,
                                GenomeCallingComponents& genome_components)
        :
        reference {genome_components.reference()},
        read_manager {genome_components.read_manager()},
        regions {genome_components.search_regions().at(contig)},
        samples {genome_components.samples()},
        caller {genome_components.caller_factory().make(contig)},
        read_buffer_size {genome_components.read_buffer_size()},
        output {output},
        progress_meter {genome_components.progress_meter()}
        {}
        
        ~ContigCallingComponents() = default;
        
        ContigCallingComponents(const ContigCallingComponents&)            = delete;
        ContigCallingComponents& operator=(const ContigCallingComponents&) = delete;
        ContigCallingComponents(ContigCallingComponents&&)                 = default;
        ContigCallingComponents& operator=(ContigCallingComponents&&)      = default;
    };
    
    bool region_has_reads(const GenomicRegion& region, ContigCallingComponents& components)
    {
        // TODO: update this to ReadManager::has_contig_reads when implemented
        return components.read_manager.get().count_reads(components.samples.get(), region) > 0;
    }
    
    auto get_call_types(const GenomeCallingComponents& components,
                        const std::vector<ContigNameType>& contigs)
    {
        CallTypeSet result {};
        
        for (const auto& contig : components.contigs_in_output_order()) {
            const auto tmp_caller = components.caller_factory().make(contig);
            
            auto caller_call_types = tmp_caller->get_call_types();
            
            result.insert(std::begin(caller_call_types), std::end(caller_call_types));
        }
        
        return result;
    }
    
    void write_caller_output_header(GenomeCallingComponents& components, const bool sites_only = false)
    {
        const auto call_types = get_call_types(components, components.contigs_in_output_order());
        
        if (sites_only) {
            components.output() << make_vcf_header({}, components.contigs_in_output_order(),
                                                   components.reference(), call_types);
        } else {
            components.output() << make_vcf_header(components.samples(),
                                                   components.contigs_in_output_order(),
                                                   components.reference(), call_types);
        }
    }
    
    void log_startup_info(const GenomeCallingComponents& components)
    {
        std::ostringstream ss {};
        
        if (components.samples().size() == 1) {
            ss << "Sample is: ";
        } else {
            ss << "Samples are: ";
        }
        
        std::transform(std::cbegin(components.samples()), std::cend(components.samples()),
                       std::ostream_iterator<std::string> {ss, " "},
                       [] (const auto& sample) -> std::string {
                           return "\"" + sample + "\"";
                       });
        
        auto str = ss.str();
        str.pop_back(); // the extra whitespace
        
        Logging::InfoLogger log {};
        
        log << str;
        
        stream(log) << "Writing calls to " << components.output().path();
    }
    
    void write_calls(std::deque<VcfRecord>&& calls, VcfWriter& out)
    {
        if (DEBUG_MODE) {
            Logging::DebugLogger log {};
            stream(log) << "Writing " << calls.size() << " calls to output";
        }
        write(calls, out);
        calls.clear();
        calls.shrink_to_fit();
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
    
    namespace
    {
        auto mapped_begin(const VcfRecord& call)
        {
            return call.pos() - 1;
        }
        
        auto mapped_end(const VcfRecord& call)
        {
            return mapped_begin(call) + static_cast<GenomicRegion::SizeType>(call.ref().size());
        }
        
        auto mapped_region(const VcfRecord& call)
        {
            return GenomicRegion {call.chrom(), mapped_begin(call), mapped_end(call)};
        }
        
        template <typename Container>
        auto encompassing_call_region(const Container& calls)
        {
            auto it = std::max_element(std::cbegin(calls), std::cend(calls),
                                       [] (const auto& lhs, const auto& rhs) {
                                           return mapped_end(lhs) < mapped_end(rhs);
                                       });
            return GenomicRegion {it->chrom(), mapped_begin(calls.front()), mapped_end(*it)};
        }
    } // namespace
    
    void buffer_connecting_calls(std::deque<VcfRecord>& calls,
                                 const GenomicRegion& next_calling_region,
                                 std::vector<VcfRecord>& buffer)
    {
        const auto it = std::find_if(std::begin(calls), std::end(calls),
                                     [&next_calling_region] (const auto& call) {
                                         return mapped_end(call) > next_calling_region.begin();
                                     });
        
        buffer.insert(std::end(buffer),
                      std::make_move_iterator(it),
                      std::make_move_iterator(std::end(calls)));
        
        calls.erase(it, std::end(calls));
    }
    
    void buffer_connecting_calls(const GenomicRegion& buffered_region,
                                 std::deque<VcfRecord>& calls,
                                 std::vector<VcfRecord>& buffer)
    {
        const auto it = std::find_if_not(std::begin(calls), std::end(calls),
                                         [&buffered_region] (const auto& call) {
                                             return mapped_begin(call) < buffered_region.end();
                                         });
        
        buffer.insert(std::end(buffer),
                      std::make_move_iterator(std::begin(calls)),
                      std::make_move_iterator(it));
        
        calls.erase(std::begin(calls), it);
    }
                
    bool is_consistent(const std::deque<VcfRecord>& merged_calls)
    {
        return true; // TODO
    }
                
    void resolve_connecting_calls(std::vector<VcfRecord>& old_connecting_calls,
                                  std::deque<VcfRecord>& calls,
                                  const ContigCallingComponents& components)
    {
        using std::begin; using std::end; using std::make_move_iterator;
        
        if (!old_connecting_calls.empty()) {
            const auto old_connecting_calls_region = encompassing_call_region(old_connecting_calls);
            
            std::vector<VcfRecord> new_connecting_calls {};
            
            buffer_connecting_calls(old_connecting_calls_region, calls, new_connecting_calls);
            
            std::deque<VcfRecord> merged_calls {};
            
            std::set_union(make_move_iterator(begin(old_connecting_calls)),
                           make_move_iterator(end(old_connecting_calls)),
                           make_move_iterator(begin(new_connecting_calls)),
                           make_move_iterator(end(new_connecting_calls)),
                           std::back_inserter(merged_calls));
            
            old_connecting_calls.clear();
            old_connecting_calls.shrink_to_fit();
            new_connecting_calls.clear();
            new_connecting_calls.shrink_to_fit();
            
            merged_calls.erase(std::unique(begin(merged_calls), end(merged_calls),
                                           [] (const auto& lhs, const auto& rhs) {
                                               return lhs.pos() == rhs.pos()
                                               && lhs.ref() == rhs.ref()
                                               && lhs.alt() == rhs.alt();
                                           }),
                               end(merged_calls));
            
            if (is_consistent(merged_calls)) {
                calls.insert(begin(calls),
                             make_move_iterator(begin(merged_calls)),
                             make_move_iterator(end(merged_calls)));
            } else {
                const auto unresolved_region = encompassing_call_region(merged_calls);
                
                merged_calls.clear();
                merged_calls.shrink_to_fit();
                
                auto new_calls = components.caller->call(unresolved_region, components.progress_meter);
                
                // TODO: we need to make sure the new calls don't contain any calls
                // outside the unresolved_region, and also possibly adjust phase regions
                // in calls past unresolved_region.
                
                calls.insert(begin(calls),
                             make_move_iterator(begin(new_calls)),
                             make_move_iterator(end(new_calls)));
            }
        }
    }
    
    void run_octopus_on_contig(ContigCallingComponents&& components)
    {
        assert(!components.regions.empty());
        
        #ifdef BENCHMARK
        init_timers();
        #endif
        
        std::deque<VcfRecord> calls;
        std::vector<VcfRecord> connecting_calls {};
        
        auto input_region = components.regions.front();
        auto subregion    = propose_call_subregion(components, input_region);
        
        auto first_input_region      = std::cbegin(components.regions);
        const auto last_input_region = std::cend(components.regions);
        
        while (!is_empty(subregion)) {
            if (DEBUG_MODE) {
                Logging::DebugLogger log {};
                stream(log) << "Processing subregion " << subregion;
            }
            
            try {
                calls = components.caller->call(subregion, components.progress_meter);
            } catch(...) {
                // TODO: which exceptions can we recover from?
                throw;
            }
            
            resolve_connecting_calls(connecting_calls, calls, components);
            
            auto next_subregion = propose_call_subregion(components, subregion, input_region);
            
            if (is_empty(next_subregion)) {
                ++first_input_region;
                if (first_input_region != last_input_region) {
                    input_region = *first_input_region;
                    next_subregion = propose_call_subregion(components, input_region);
                }
            }
            
            assert(connecting_calls.empty());
            
            buffer_connecting_calls(calls, next_subregion, connecting_calls);
            
            try {
                write_calls(std::move(calls), components.output);
            } catch(...) {
                // TODO: which exceptions can we recover from?
                throw;
            }
            
            subregion = std::move(next_subregion);
        }
        
        #ifdef BENCHMARK
        print_all_timers();
        #endif
    }
    
    void run_octopus_single_threaded(GenomeCallingComponents& components)
    {
        components.progress_meter().start();
        
        for (const auto& contig : components.contigs_in_output_order()) {
            run_octopus_on_contig(ContigCallingComponents {contig, components});
        }
        
        components.progress_meter().stop();
    }
    
    VcfWriter create_unique_temp_output_file(const GenomicRegion& region,
                                             const GenomeCallingComponents& components)
    {
        auto path = *components.temp_directory();
        
        const auto& contig = region.contig_name();
        const auto begin   = std::to_string(region.begin());
        const auto end     = std::to_string(region.end());
        
        auto file_name = fs::path {contig + "_" + begin + "-" + end + "_temp.bcf"};
        
        path /= file_name;
        
        const auto call_types = get_call_types(components, {region.contig_name()});
        
        auto header = make_vcf_header(components.samples(), contig, components.reference(), call_types);
        
        return VcfWriter {std::move(path), std::move(header)};
    }
    
    VcfWriter create_unique_temp_output_file(const GenomicRegion::ContigNameType& contig,
                                             const GenomeCallingComponents& components)
    {
        return create_unique_temp_output_file(components.reference().contig_region(contig),
                                              components);
    }
    
    using TempVcfWriterMap = std::unordered_map<ContigNameType, VcfWriter>;
    
    TempVcfWriterMap make_temp_writers(const GenomeCallingComponents& components)
    {
        if (!components.temp_directory()) {
            throw std::runtime_error {"Could not make temp writers"};
        }
        
        TempVcfWriterMap result {};
        result.reserve(components.contigs_in_output_order().size());
        
        for (const auto& contig : components.contigs_in_output_order()) {
            result.emplace(contig, create_unique_temp_output_file(contig, components));
        }
        
        return result;
    }
    
    struct Task : public Mappable<Task>
    {
        enum class ExecutionPolicy { Threaded, VectorThreaded, None };
        
        Task() = delete;
        
        Task(GenomicRegion region, ExecutionPolicy policy = ExecutionPolicy::None)
        : region {std::move(region)}, policy {policy} {};
        
        GenomicRegion region;
        ExecutionPolicy policy;
        
        const GenomicRegion& mapped_region() const noexcept { return region; }
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
    using ContigTaskMap = std::map<ContigNameType, TaskQueue, ContigOrder>;
    
    TaskQueue divide_work_into_tasks(const ContigCallingComponents& components,
                                     const Task::ExecutionPolicy policy)
    {
        TaskQueue result {};
        
        if (components.regions.empty()) return result;
        
        for (const auto& region : components.regions) {
            auto subregion = propose_call_subregion(components, region);
            
            while (!is_empty(subregion)) {
                result.emplace(subregion, policy);
                subregion = propose_call_subregion(components, subregion, region);
            }
            
            if (result.empty()) {
                result.emplace(region, policy);
            }
        }
        
        return result;
    }
    
    Task::ExecutionPolicy make_execution_policy(const GenomeCallingComponents& components)
    {
        if (components.num_threads()) {
            return Task::ExecutionPolicy::None;
        }
        return Task::ExecutionPolicy::Threaded;
    }
    
    ContigTaskMap make_tasks(GenomeCallingComponents& components, const unsigned num_threads)
    {
        const auto policy = make_execution_policy(components);
        
        ContigTaskMap result {ContigOrder {components.contigs_in_output_order()}};
        
        for (const auto& contig : components.contigs_in_output_order()) {
            ContigCallingComponents contig_components {contig, components};
            contig_components.read_buffer_size /= num_threads;
            result.emplace(contig, divide_work_into_tasks(contig_components, policy));
        }
        
        return result;
    }
    
    unsigned calculate_num_task_threads(const GenomeCallingComponents& components)
    {
        if (components.num_threads()) {
            return *components.num_threads();
        }
        
        // TODO: come up with a better calculation
        
        const auto num_hardware_threads = std::thread::hardware_concurrency();
        
        if (num_hardware_threads > 0) return num_hardware_threads;
        
        Logging::WarningLogger log {};
        log << "Unable to detect the number of system threads,"
            " it may be better to run with a user number if the number of cores is known";
        
        const auto num_files = components.read_manager().num_files();
        
        return std::min(num_files, 8u);
    }
    
    Task pop(ContigTaskMap& tasks)
    {
        assert(!tasks.empty());
        auto it = std::begin(tasks);
        assert(!it->second.empty());
        const auto result = it->second.front();
        it->second.pop();
        if (it->second.empty()) {
            tasks.erase(it);
        }
        return result;
    }
    
    template<typename R>
    bool is_ready(const std::future<R>& f)
    {
        return f.valid() && f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
    }
    
    struct CompletedTask : public Task, public Mappable<CompletedTask>
    {
        CompletedTask(Task task, std::deque<VcfRecord>&& calls)
        : Task {std::move(task)}, calls {std::move(calls)} {}
        
        std::deque<VcfRecord> calls;
    };
    
    struct SyncPacket
    {
        std::condition_variable cv;
        std::mutex mutex;
        std::atomic<unsigned> count_finished;
    };
    
    auto run(Task task, GenomeCallingComponents& components, SyncPacket& sync)
    {
        if (DEBUG_MODE) {
            Logging::DebugLogger log {};
            stream(log) << "Spawning task with region " << task.region;
        }
        
        ContigCallingComponents contig_components {contig_name(task), components};
        
        return std::async(std::launch::async,
                          [task = std::move(task), components = std::move(contig_components),
                           &sync] () -> CompletedTask {
                              CompletedTask result {
                                  task,
                                  components.caller->call(task.region, components.progress_meter)
                              };
                              
                              std::unique_lock<std::mutex> lk {sync.mutex};
                              ++sync.count_finished;
                              lk.unlock();
                              sync.cv.notify_all();
                              
                              return result;
                          });
    }
    
    template <typename K>
    auto extract(std::unordered_map<K, VcfWriter>& writers)
    {
        std::vector<VcfWriter> result {};
        result.reserve(writers.size());
        
        for (auto& p : writers) result.push_back(std::move(p.second));
        
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
        
        // TODO: start running tasks as soon as they are added into the task map,
        // this will require some additional syncronisation on the task map
        
        auto tasks = make_tasks(components, num_task_threads);
        
        TempVcfWriterMap temp_writers;
        try {
            temp_writers = make_temp_writers(components);
        } catch (...) {
            Logging::ErrorLogger log {};
            log << "Could not make required temporary files";
            return;
        }
        
        std::vector<std::future<CompletedTask>> futures(num_task_threads);
        
        ContigTaskMap pending_tasks {ContigOrder {components.contigs_in_output_order()}};
        
        MappableSetMap<ContigNameType, CompletedTask> buffered_tasks {};
        
        SyncPacket sync {};
        
        components.progress_meter().start();
        
        while (!tasks.empty()) {
            for (auto& future : futures) {
                if (is_ready(future)) {
                    try {
                        std::lock_guard<std::mutex> lk {sync.mutex};
                        
                        auto result = future.get();
                        
                        const auto& contig = contig_name(result);
                        
                        assert(pending_tasks.count(contig) == 1);
                        
                        auto& contig_pending_tasks = pending_tasks.at(contig);
                        
                        if (is_same_region(result, contig_pending_tasks.front())) {
                            if (DEBUG_MODE) {
                                Logging::DebugLogger log {};
                                stream(log) << "Writing completed task " << result.region;
                            }
                            
                            write_calls(std::move(result.calls), temp_writers.at(contig_name(result)));
                            
                            contig_pending_tasks.pop();
                            
                            if (buffered_tasks.count(contig) == 1) {
                                auto& contig_buffered_tasks = buffered_tasks.at(contig);
                                
                                while (!contig_pending_tasks.empty()) {
                                    if (contig_buffered_tasks.has_contained(contig_pending_tasks.front())) {
                                        if (DEBUG_MODE) {
                                            Logging::DebugLogger log {};
                                            stream(log) << "Writing completed buffered task " << contig_pending_tasks.front().region;
                                        }
                                        
                                        auto buffered_task = contig_buffered_tasks.contained_range(contig_pending_tasks.front()).front();
                                        
                                        contig_buffered_tasks.erase(buffered_task);
                                        
                                        write_calls(std::move(buffered_task.calls),
                                                    temp_writers.at(contig_name(result)));
                                        
                                        contig_pending_tasks.pop();
                                    } else {
                                        break;
                                    }
                                }
                            }
                        } else {
                            if (DEBUG_MODE) {
                                Logging::DebugLogger log {};
                                stream(log) << "Buffering completed task " << result.region;
                            }
                            
                            buffered_tasks[contig].insert(std::move(result));
                        }
                        
                        --sync.count_finished;
                    } catch (...) {
                        // TODO: which exceptions can we recover from?
                        throw;
                    }
                }
                
                if (!tasks.empty() && !future.valid()) {
                    std::lock_guard<std::mutex> lk {sync.mutex};
                    auto task = pop(tasks);
                    
                    future = run(task, components, sync);
                    
                    pending_tasks[task.region.contig_name()].push(std::move(task));
                }
            }
            
            if (sync.count_finished == 0) {
                std::unique_lock<std::mutex> lk {sync.mutex};
                while (sync.count_finished == 0) { // for spurious wakeup
                    sync.cv.wait(lk);
                }
            }
        }
        
        pending_tasks.clear();
        
        futures.erase(std::remove_if(std::begin(futures), std::end(futures),
                                     [] (const auto& f) { return !f.valid(); }),
                      std::end(futures));
        
        std::deque<CompletedTask> remaining_tasks {};
        
        std::transform(std::begin(futures), std::end(futures), std::back_inserter(remaining_tasks),
                       [] (auto& fut) {
                           try {
                               return fut.get();
                           } catch (...) {
                               // TODO: which exceptions can we recover from?
                               throw;
                           }
                       });
        
        futures.clear();
        futures.shrink_to_fit();
        
        for (auto& p : buffered_tasks) {
            remaining_tasks.insert(std::end(remaining_tasks),
                                   std::make_move_iterator(std::begin(p.second)),
                                   std::make_move_iterator(std::end(p.second)));
            p.second.clear();
            p.second.shrink_to_fit();
        }
        
        buffered_tasks.clear();
        
        std::sort(std::begin(remaining_tasks), std::end(remaining_tasks),
                  [] (const auto& lhs, const auto& rhs) {
                      if (contig_name(lhs) == contig_name(rhs)) {
                          return lhs < rhs;
                      }
                      return contig_name(lhs) < contig_name(rhs);
                  });
        
        for (auto& task : remaining_tasks) {
            if (DEBUG_MODE) {
                Logging::DebugLogger log {};
                stream(log) << "Writing remaining completed task " << task.region;
            }
            write_calls(std::move(task.calls), temp_writers.at(contig_name(task)));
        }
        
        remaining_tasks.clear();
        remaining_tasks.shrink_to_fit();
        
        components.progress_meter().stop();
        
        auto results = convert_temp_writers(temp_writers);
        
        index_vcfs(results);
        
        merge(results, components.output(), components.contigs_in_output_order());
    }
    } // namespace
    
    bool is_multithreaded(const GenomeCallingComponents& components)
    {
        return !components.num_threads() || *components.num_threads() > 1;
    }
    
    auto make_filter_read_pipe(const GenomeCallingComponents& components)
    {
        using std::make_unique;
        
        ReadTransform transform {};
        transform.register_transform(ReadTransforms::MaskSoftClipped {});
        
        ReadFilterer filter {};
        filter.register_filter(make_unique<ReadFilters::HasValidQualities>());
        filter.register_filter(make_unique<ReadFilters::HasWellFormedCigar>());
        filter.register_filter(make_unique<ReadFilters::IsMapped>());
        filter.register_filter(make_unique<ReadFilters::IsNotMarkedQcFail>());
        filter.register_filter(make_unique<ReadFilters::IsNotMarkedDuplicate>());
        filter.register_filter(make_unique<ReadFilters::IsNotDuplicate<ReadFilterer::BidirIt>>());
        filter.register_filter(make_unique<ReadFilters::IsProperTemplate>());
        
        return ReadPipe {
            components.read_manager(), std::move(transform), std::move(filter), boost::none,
            components.samples()
        };
    }
    
    auto get_filtered_path(const fs::path& unfiltered_output_path)
    {
        const std::string identifier {"filtered"};
        const fs::path new_stem {
            unfiltered_output_path.stem().string()
            + "_" + identifier
            + unfiltered_output_path.extension().string()
        };
        return unfiltered_output_path.parent_path() / new_stem;
    }
                
    auto get_filtered_path(const GenomeCallingComponents& components)
    {
        return get_filtered_path(components.output().path());
    }
    
    void filter_calls(const GenomeCallingComponents& components)
    {
        //assert(!components.output().is_open());
        
        const VcfReader calls {components.output().path()};
        //const VcfReader calls {"/Users/danielcooke/Genomics/octopus_test/octopus_calls4.vcf"};
        
        const auto filtered_path = get_filtered_path(components);
        //const auto filtered_path = get_filtered_path(calls.path());
        
        VcfWriter filtered_calls {filtered_path};
        
        const auto read_pipe = make_filter_read_pipe(components);
        
        const VariantCallFilter filter {components.reference(), read_pipe};
        
        VariantCallFilter::RegionMap regions {
            [contigs = components.contigs_in_output_order()]
            (const auto& lhs, const auto& rhs) -> bool {
                std::array<std::reference_wrapper<const std::decay_t<decltype(lhs)>>, 2> values {lhs, rhs};
                return *std::find_first_of(std::cbegin(contigs), std::cend(contigs),
                                           std::cbegin(values), std::cend(values),
                                           [] (const auto& lhs, const auto& rhs) {
                                               return lhs == rhs.get();
                                           }) == lhs;
            }
        };
        
        for (const auto& p : components.search_regions()) {
            regions.emplace(std::piecewise_construct,
                            std::forward_as_tuple(p.first),
                            std::forward_as_tuple(std::initializer_list<GenomicRegion> {encompassing_region(p.second)}));
                            //std::forward_as_tuple(std::cbegin(p.second), std::cend(p.second)));
        }
        
        filter.filter(calls, filtered_calls, regions);
    }
                
    void convert_to_legacy(const VcfReader& src, VcfWriter& dst)
    {
        if (!dst.is_header_written()) {
            dst << src.fetch_header();
        }
        
        auto calls = src.fetch_records();
        
        const auto samples = src.fetch_header().samples();
        
        const static char Deleted {'*'};
        const static std::string Missing {"."};
        
        const auto has_deleted = [] (const auto& allele) {
            return std::find(std::cbegin(allele), std::cend(allele), Deleted) != std::cend(allele);
        };
        
        for (const auto& call : calls) {
            const auto& alt = call.alt();
            
            VcfRecord::Builder cb {call};
            
            const auto it = std::find_if(std::cbegin(alt), std::cend(alt), has_deleted);
            
            if (it != std::cend(alt)) {
                const auto i = std::distance(std::cbegin(alt), it);
                
                auto new_alt = alt;
                
                new_alt.erase(std::next(std::begin(new_alt), i));
                
                cb.set_alt(std::move(new_alt));
            }
            
            for (const auto& sample : samples) {
                const auto& gt = call.get_sample_value(sample, "GT");
                
                const auto it2 = std::find_if(std::cbegin(gt), std::cend(gt), has_deleted);
                const auto it3 = std::find(std::cbegin(gt), std::cend(gt), Missing);
                
                const auto& ref = call.ref();
                
                if (it2 != std::cend(gt) || it3 != std::cend(gt)) {
                    auto new_gt = gt;
                    
                    std::replace_if(std::begin(new_gt), std::end(new_gt), has_deleted, ref);
                    
                    std::replace(std::begin(new_gt), std::end(new_gt), Missing, ref);
                    
                    auto phasing = VcfRecord::Builder::Phasing::Phased;
                    
                    if (!call.is_sample_phased(sample)) {
                        phasing = VcfRecord::Builder::Phasing::Unphased;
                    }
                    
                    cb.set_genotype(sample, std::move(new_gt), phasing);
                }
            }
            
            dst << cb.build_once();
        }
    }
    
    auto get_legacy_path(const fs::path& base)
    {
        const std::string identifier {"legacy"};
        const fs::path new_stem {base.stem().string() + "." + identifier + base.extension().string()};
        return base.parent_path() / new_stem;
    }
                
    auto get_legacy_path(const GenomeCallingComponents& components)
    {
        return get_legacy_path(components.output().path());
    }
    
    void run_octopus(po::variables_map& options)
    {
        DEBUG_MODE = Options::is_debug_mode(options);
        TRACE_MODE = Options::is_trace_mode(options);
        
        log_startup();
        
        Logging::InfoLogger log {};
        
        const auto start = std::chrono::system_clock::now();
        
        auto end = start;
        
        unsigned search_size {0};
        
        // open scope to ensure calling components are destroyed before end message
        {
            auto components = collate_genome_calling_components(options);
            
            if (!components) return;
            
            end = std::chrono::system_clock::now();
            
            stream(log) << "Done initialising calling components in " << TimeInterval {start, end};
            
            log_startup_info(*components);
            
            write_caller_output_header(*components, Options::call_sites_only(options));
            
            options.clear();
            
            try {
                if (is_multithreaded(*components)) {
                    if (DEBUG_MODE) {
                        Logging::WarningLogger warn_log {};
                        warn_log << "Running in parallel mode can make debug log difficult to interpret";
                    }
                    
                    run_octopus_multi_threaded(*components);
                } else {
                    run_octopus_single_threaded(*components);
                }
                
                components->output().close();
            } catch (const std::exception& e) {
                Logging::FatalLogger lg {};
                stream(lg) << "Encountered exception '" << e.what() << "'. Attempting to cleanup...";
                
                cleanup(*components);
                
                if (DEBUG_MODE) {
                    stream(lg) << "Cleanup successful. Please send log file to dcooke@well.ox.ac.uk";
                } else {
                    stream(lg) << "Cleanup successful. Please re-run in debug mode (option --debug) and send"
                                    " log file to " << Octopus_bug_email;
                }
                
                return;
            }
            
            filter_calls(*components);
            
            const auto filtered_path = get_filtered_path(*components);
            
            const VcfReader vcf {filtered_path};
            VcfWriter out {get_legacy_path(filtered_path)};
//            const VcfReader vcf {"/Users/danielcooke/Genomics/octopus_test/octopus_calls4_filtered.vcf"};
//            VcfWriter out {"/Users/danielcooke/Genomics/octopus_test/octopus_calls4_filtered.legacy.vcf.gz"};
            convert_to_legacy(vcf, out);
            
            cleanup(*components);
            
            end = std::chrono::system_clock::now();
            
            search_size= calculate_total_search_size(components->search_regions());
        }
        
        stream(log) << "Finished processing " << search_size << "bp, total runtime " << TimeInterval {start, end};
    }
    
} // namespace Octopus


