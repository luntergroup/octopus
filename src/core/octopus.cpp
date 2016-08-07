// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

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
#include <set>
#include <typeinfo>
#include <cassert>

#include <boost/program_options.hpp>
#include <boost/optional.hpp>

#include <config/common.hpp>
#include <config/option_collation.hpp>
#include <basics/genomic_region.hpp>
#include <concepts/mappable.hpp>
#include <containers/mappable_flat_multi_set.hpp>
#include <containers/mappable_map.hpp>
#include <io/reference/reference_genome.hpp>
#include <io/read/read_manager.hpp>
#include <readpipe/read_pipe_fwd.hpp>
#include <utils/read_stats.hpp>
#include <config/octopus_vcf.hpp>
#include <core/callers/caller_factory.hpp>
#include <core/callers/caller.hpp>
#include <utils/maths.hpp>
#include <logging/progress_meter.hpp>
#include <logging/logging.hpp>
#include <core/callers/utils/vcf_header_factory.hpp>
#include <io/variant/vcf.hpp>
#include <utils/timing.hpp>

#include "timers.hpp" // BENCHMARK

namespace octopus {

using options::OptionMap;
using logging::get_debug_log;

namespace {

template <typename T>
std::size_t index_of(const std::vector<T>& elements, const T& value)
{
    return std::distance(std::cbegin(elements), std::find(std::cbegin(elements), std::cend(elements), value));
}

auto get_contigs(const InputRegionMap& regions, const ReferenceGenome& reference,
                 const options::ContigOutputOrder order)
{
    using options::ContigOutputOrder;
    
    using ContigName = GenomicRegion::ContigName;
    
    std::vector<ContigName> result {};
    result.reserve(regions.size());
    
    std::transform(std::cbegin(regions), std::cend(regions), std::back_inserter(result),
                   [] (const auto& p) { return p.first; });
    
    switch (order) {
        case ContigOutputOrder::LexicographicalAscending:
            std::sort(std::begin(result), std::end(result));
            break;
        case ContigOutputOrder::LexicographicalDescending:
            std::sort(std::begin(result), std::end(result), std::greater<ContigName>());
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

template <typename S>
void print_input_regions(S&& stream, const InputRegionMap& regions)
{
    stream << "All input regions:" << '\n';
    for (const auto& p : regions) {
        stream << "Contig " << p.first << '\n';
        for (const auto& region : p.second) {
            stream << region << ' ';
        }
        stream << '\n';
    }
}

void print_input_regions(const InputRegionMap& regions)
{
    print_input_regions(std::cout, regions);
}

template <typename Container>
bool is_in_file_samples(const SampleName& sample, const Container& file_samples)
{
    return std::find(std::cbegin(file_samples), std::cend(file_samples), sample)
                != std::cend(file_samples);
}

std::vector<SampleName> extract_samples(const OptionMap& options,
                                          const ReadManager& read_manager)
{
    auto user_samples = options::get_user_samples(options);
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
            std::transform(it, std::end(*user_samples), std::ostream_iterator<SampleName>(ss, ", "),
                           [] (auto sample) { return "'" + sample + "'"; });
            if (num_not_found == 1) {
                ss << "is";
            } else {
                ss << "are";
            }
            ss << " not present in any of the read files";
            logging::WarningLogger log {};
            log << ss.str();
            user_samples->erase(it, std::end(*user_samples));
        }
        
        return *user_samples;
    } else {
        return file_samples;
    }
}

using CallTypeSet = std::set<std::type_index>;

VcfHeader make_vcf_header(const std::vector<SampleName>& samples,
                          const std::vector<GenomicRegion::ContigName>& contigs,
                          const ReferenceGenome& reference,
                          const CallTypeSet& call_types)
{
    auto builder = Vcf::make_octopus_header_template().set_samples(samples);
    
    for (const auto& contig : contigs) {
        builder.add_contig(contig, {{"length", std::to_string(reference.contig_size(contig))}});
    }
    
    builder.add_basic_field("reference", reference.name());
    builder.add_structured_field("octopus", {{"some", "option"}});
    
    VcfHeaderFactory factory {};
    
    for (const auto& type : call_types) {
        factory.register_call_type(type);
    }
    
    factory.annotate(builder);
    
    return builder.build_once();
}

VcfHeader make_vcf_header(const std::vector<SampleName>& samples,
                          const GenomicRegion::ContigName& contig,
                          const ReferenceGenome& reference,
                          const CallTypeSet& call_types)
{
    return make_vcf_header(samples, std::vector<GenomicRegion::ContigName> {contig},
                           reference, call_types);
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
        + sequence_size(read) * sizeof(AlignedRead::BaseQuality)
        + read.cigar().size() * sizeof(CigarOperation)
        + contig_name(read).size()
        + (read.has_other_segment() ? sizeof(AlignedRead::Segment) : 0);
}

auto estimate_mean_read_size(const std::vector<SampleName>& samples,
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
        logging::WarningLogger log {};
        log << "Could not estimate read size from data, resorting to default";
        
        return sizeof(AlignedRead) + 300;
    }
    
    return static_cast<std::size_t>(maths::mean(read_sizes) + maths::stdev(read_sizes));
}

std::size_t calculate_max_num_reads(const std::size_t max_buffer_bytes,
                                    const std::vector<SampleName>& samples,
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
    using Samples = std::vector<SampleName>;
    using Contigs = std::vector<GenomicRegion::ContigName>;
    
    GenomeCallingComponents() = delete;
    
    GenomeCallingComponents(ReferenceGenome&& reference, ReadManager&& read_manager,
                            VcfWriter&& output, const OptionMap& options)
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
    
    const CallerFactory& caller_factory() const noexcept
    {
        return components_.caller_factory;
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
                   VcfWriter&& output, const OptionMap& options)
        :
        reference {std::move(reference)},
        read_manager {std::move(read_manager)},
        samples {extract_samples(options, this->read_manager)},
        regions {options::get_search_regions(options, this->reference)},
        contigs_in_output_order {get_contigs(this->regions, this->reference,
                                             options::get_contig_output_order(options))},
        read_pipe {options::make_read_pipe(this->read_manager, this->samples, options)},
        caller_factory {options::make_caller_factory(this->reference,
                                                                     this->read_pipe,
                                                                     this->regions,
                                                                     options)},
        output {std::move(output)},
        num_threads {options::get_num_threads(options)},
        read_buffer_size {},
        temp_directory {((!num_threads || *num_threads > 1) ? options::create_temp_file_directory(options) : boost::none)},
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
                read_buffer_size = calculate_max_num_reads(options::get_target_read_buffer_size(options),
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
        CallerFactory caller_factory;
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
        components_.caller_factory.set_reference(components_.reference);
        components_.caller_factory.set_read_pipe(components_.read_pipe);
    }
};

bool are_components_valid(const GenomeCallingComponents& components)
{
    logging::FatalLogger log {};
    
    if (components.samples().empty()) {
        log << "No samples detected - at least one is required for calling";
        return false;
    }
    
    if (components.search_regions().empty()) {
        log << "There are no input regions - at least one is required for calling";
        return false;
    }
    
    return true;
}

void cleanup(GenomeCallingComponents& components) noexcept
{
    logging::InfoLogger log {};
    if (components.temp_directory()) {
        try {
            const auto num_files_removed = fs::remove_all(*components.temp_directory());
            stream(log) << "Removed " << num_files_removed << " temporary files";
        } catch (const std::exception& e) {
            stream(log) << "Cleanup failed with exception: " << e.what();
        }
    }
}

boost::optional<GenomeCallingComponents> collate_genome_calling_components(const OptionMap& options)
{
    auto reference = options::make_reference(options);
    
    auto read_manager = options::make_read_manager(options);
    
    auto output = options::make_output_vcf_writer(options);
    
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
}

struct ContigCallingComponents
{
    std::reference_wrapper<const ReferenceGenome> reference;
    std::reference_wrapper<ReadManager> read_manager;
    const InputRegionMap::mapped_type regions;
    std::reference_wrapper<const std::vector<SampleName>> samples;
    std::unique_ptr<const Caller> caller;
    std::size_t read_buffer_size;
    std::reference_wrapper<VcfWriter> output;
    std::reference_wrapper<ProgressMeter> progress_meter;
    
    ContigCallingComponents() = delete;
    
    ContigCallingComponents(const GenomicRegion::ContigName& contig,
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
    
    ContigCallingComponents(const GenomicRegion::ContigName& contig, VcfWriter& output,
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

bool has_reads(const GenomicRegion& region, ContigCallingComponents& components)
{
    return components.read_manager.get().has_reads(components.samples.get(), region);
}

auto get_call_types(const GenomeCallingComponents& components,
                    const std::vector<ContigName>& contigs)
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
    
    logging::InfoLogger log {};
    
    log << str;
    
    stream(log) << "Writing calls to " << components.output().path();
}

void write_calls(std::deque<VcfRecord>&& calls, VcfWriter& out)
{
    static auto debug_log = get_debug_log();
    if (debug_log) stream(*debug_log) << "Writing " << calls.size() << " calls to output";
    write(calls, out);
    calls.clear();
    calls.shrink_to_fit();
}

auto find_max_window(const ContigCallingComponents& components,
                     const GenomicRegion& remaining_call_region)
{
    const auto& rm = components.read_manager.get();
    
    if (!rm.has_reads(components.samples.get(), remaining_call_region)) {
        return remaining_call_region;
    }
    
    auto result = rm.find_covered_subregion(components.samples, remaining_call_region,
                                            components.read_buffer_size);
    
    if (ends_before(result, remaining_call_region)) {
        auto rest = right_overhang_region(remaining_call_region, result);
        
        if (!rm.has_reads(components.samples.get(), rest)) {
            result = remaining_call_region;
        }
    }
    
    return result;
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
        return mapped_begin(call) + static_cast<GenomicRegion::Position>(call.ref().size());
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
    // TODO: refactor to use connection resolution developed for multithreaded version
    
    static auto debug_log = get_debug_log();
    
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
    
    while (first_input_region != last_input_region && !is_empty(subregion)) {
        if (debug_log) stream(*debug_log) << "Processing subregion " << subregion;
        
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

VcfWriter create_unique_temp_output_file(const GenomicRegion::ContigName& contig,
                                         const GenomeCallingComponents& components)
{
    return create_unique_temp_output_file(components.reference().contig_region(contig),
                                          components);
}

using TempVcfWriterMap = std::unordered_map<ContigName, VcfWriter>;

TempVcfWriterMap make_temp_vcf_writers(const GenomeCallingComponents& components)
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

std::ostream& operator<<(std::ostream& os, const Task& task)
{
    os << task.region;
    return os;
}

struct ContigOrder
{
    using ContigName = GenomicRegion::ContigName;
    
    template <typename Container>
    ContigOrder(const Container& contigs)
    : contigs_ {std::cbegin(contigs), std::cend(contigs)} {}
    
    bool operator()(const ContigName& lhs, const ContigName& rhs) const
    {
        const auto it1 = std::find(std::cbegin(contigs_), std::cend(contigs_), lhs);
        const auto it2 = std::find(std::cbegin(contigs_), std::cend(contigs_), rhs);
        return it1 < it2;
    }
    
private:
    std::vector<ContigName> contigs_;
};

using TaskQueue = std::queue<Task>;
using TaskMap   = std::map<ContigName, TaskQueue, ContigOrder>;

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

TaskMap make_tasks(GenomeCallingComponents& components, const unsigned num_threads)
{
    const auto policy = make_execution_policy(components);
    
    TaskMap result {ContigOrder {components.contigs_in_output_order()}};
    
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
    
    logging::WarningLogger log {};
    log << "Unable to detect the number of system threads,"
        " it may be better to run with a user number if the number of cores is known";
    
    const auto num_files = components.read_manager().num_files();
    
    return std::min(num_files, 8u);
}

Task pop(TaskMap& tasks)
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
    CompletedTask(Task task) : Task {std::move(task)} {}
    CompletedTask(Task task, std::deque<VcfRecord>&& calls)
    : Task {std::move(task)}, calls {std::move(calls)} {}
    
    std::deque<VcfRecord> calls;
};

std::ostream& operator<<(std::ostream& os, const CompletedTask& task)
{
    os << task.region;
    return os;
}

struct SyncPacket
{
    std::condition_variable cv;
    std::mutex mutex;
    std::atomic<unsigned> count_finished;
};

auto run(Task task, ContigCallingComponents components, SyncPacket& sync)
{
    static auto debug_log = get_debug_log();
    
    if (debug_log) stream(*debug_log) << "Spawning task with region " << task.region;
    
    return std::async(std::launch::async,
                      [task = std::move(task), components = std::move(components),
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

using CompletedTaskMap = std::map<ContigName, std::map<ContigRegion, CompletedTask>>;

using HoldbackTask = boost::optional<std::reference_wrapper<const CompletedTask>>;

auto get_writable_completed_tasks(CompletedTask&& task, CompletedTaskMap::mapped_type& buffered_tasks,
                                  TaskQueue& running_tasks, HoldbackTask& holdback)
{
    std::deque<CompletedTask> result {std::move(task)};
    
    while (!running_tasks.empty()) {
        const auto it = buffered_tasks.find(contig_region(running_tasks.front()));
        
        if (it != std::end(buffered_tasks)) {
            result.push_back(std::move(it->second));
            buffered_tasks.erase(it);
            running_tasks.pop();
        } else {
            break;
        }
    }
    
    if (holdback) {
        const auto it = buffered_tasks.find(contig_region(holdback->get()));
        assert(it->second < result.front());
        result.push_front(std::move(it->second));
        buffered_tasks.erase(it);
        holdback = boost::none;
    }
    
    return result;
}

auto encompassing_call_region(const CompletedTask& task)
{
    return encompassing_call_region(task.calls);
}

using ContigCallingComponentFactory    = std::function<ContigCallingComponents()>;
using ContigCallingComponentFactoryMap = std::map<ContigName, ContigCallingComponentFactory>;

auto make_contig_calling_component_factory_map(GenomeCallingComponents& components)
{
    ContigCallingComponentFactoryMap result {};
    
    for (const auto& contig : components.contigs_in_output_order()) {
        result.emplace(contig, [&components, contig] () -> ContigCallingComponents
                       { return ContigCallingComponents {contig, components}; });
    }
    
    return result;
}

auto find_first_lhs_connecting(CompletedTask& lhs, const GenomicRegion& rhs_region)
{
    const auto rhs_begin = mapped_begin(rhs_region);
    return std::find_if(std::begin(lhs.calls), std::end(lhs.calls),
                        [&rhs_begin] (const auto& call) {
                            return mapped_end(call) > rhs_begin;
                        });
}

auto find_last_rhs_connecting(const GenomicRegion& lhs_region, CompletedTask& rhs)
{
    const auto lhs_end = mapped_end(lhs_region);
    return std::find_if_not(std::begin(rhs.calls), std::end(rhs.calls),
                            [&lhs_end] (const auto& call) {
                                return mapped_begin(call) < lhs_end;
                            });
}

void resolve_connecting_calls(CompletedTask& lhs, CompletedTask& rhs,
                              const ContigCallingComponentFactory& calling_components)
{
    auto first_lhs_connecting = find_first_lhs_connecting(lhs, encompassing_call_region(rhs));
    auto last_rhs_connecting  = find_last_rhs_connecting(encompassing_call_region(lhs), rhs);
    
    using std::begin; using std::end; using std::make_move_iterator;
    
    if (first_lhs_connecting == end(lhs.calls) && last_rhs_connecting == begin(rhs.calls)) {
        return;
    }
    
    static auto debug_log = get_debug_log();
    
    if (debug_log) stream(*debug_log) << "Resolving connecting calls between tasks " << lhs << " & " << rhs;
    
    std::deque<VcfRecord> merged_calls {};
    
    std::set_union(make_move_iterator(first_lhs_connecting),
                   make_move_iterator(end(lhs.calls)),
                   make_move_iterator(begin(rhs.calls)),
                   make_move_iterator(last_rhs_connecting),
                   std::back_inserter(merged_calls));
    
    lhs.calls.erase(first_lhs_connecting, end(lhs.calls));
    rhs.calls.erase(begin(rhs.calls), last_rhs_connecting);
    
    merged_calls.erase(std::unique(begin(merged_calls), end(merged_calls),
                                   [] (const auto& lhs, const auto& rhs) {
                                       return lhs.pos() == rhs.pos()
                                       && lhs.ref() == rhs.ref()
                                       && lhs.alt() == rhs.alt();
                                   }),
                       end(merged_calls));
    
    if (is_consistent(merged_calls)) {
        rhs.calls.insert(begin(rhs.calls),
                         make_move_iterator(begin(merged_calls)),
                         make_move_iterator(end(merged_calls)));
    } else {
        const auto unresolved_region = encompassing_region(mapped_region(merged_calls.front()),
                                                           mapped_region(merged_calls.back()));
        
        const auto components = calling_components();
        
        auto num_unresolved_region_reads = components.read_manager.get().count_reads(components.samples, unresolved_region);
        
        if (num_unresolved_region_reads <= components.read_buffer_size) {
            merged_calls.clear();
            merged_calls.shrink_to_fit();
            
            if (debug_log) {
                stream(*debug_log) << "Calls are inconsistent in connecting region " << unresolved_region
                << ". Recalling the region";
            }
            
            auto resolved_calls = components.caller->call(unresolved_region, components.progress_meter);
            
            if (!resolved_calls.empty()) {
                if (!contains(unresolved_region, encompassing_call_region(resolved_calls))) {
                    // TODO
                }
                
                // TODO: we may need to adjust phase regions in calls past the unresolved_region
                
                rhs.calls.insert(begin(rhs.calls),
                                 make_move_iterator(begin(resolved_calls)),
                                 make_move_iterator(end(resolved_calls)));
            }
        } else {
            // TODO: we could try to manually resolve the calls. Very difficult.
            logging::WarningLogger log {};
            stream(log) << "Skipping region " << unresolved_region << " as there are too many reads"
                        " to analyse the whole region, and partitions give inconsistent calls";
        }
    }
}

void resolve_connecting_calls(std::deque<CompletedTask>& adjacent_tasks,
                              const ContigCallingComponentFactory& calling_components)
{
    assert(std::is_sorted(std::cbegin(adjacent_tasks), std::cend(adjacent_tasks)));
    auto lhs = std::begin(adjacent_tasks);
    std::for_each(std::next(lhs), std::end(adjacent_tasks),
                  [&] (auto& rhs) { resolve_connecting_calls(*lhs++, rhs, calling_components); });
}

void write(std::deque<CompletedTask>&& tasks, VcfWriter& temp_vcf)
{
    static auto debug_log = get_debug_log();
    for (auto&& task : tasks) {
        if (debug_log) stream(*debug_log) << "Writing completed task " << task;
        write_calls(std::move(task.calls), temp_vcf);
    }
}

void write_or_buffer(CompletedTask&& task, CompletedTaskMap::mapped_type& buffered_tasks,
                     TaskQueue& running_tasks, HoldbackTask& holdback,
                     VcfWriter& temp_vcf, const ContigCallingComponentFactory& calling_components)
{
    static auto debug_log = get_debug_log();
    
    if (is_same_region(task, running_tasks.front())) {
        running_tasks.pop();
        
        auto writable_tasks = get_writable_completed_tasks(std::move(task), buffered_tasks,
                                                           running_tasks, holdback);
        
        assert(holdback == boost::none);
        
        resolve_connecting_calls(writable_tasks, calling_components);
        
        // keep the last task buffered to enable connection resolution when the next task finishes
        
        const auto p = buffered_tasks.emplace(contig_region(writable_tasks.back()),
                                              std::move(writable_tasks.back()));
        
        assert(p.second);
        holdback = p.first->second;
        
        if (debug_log) stream(*debug_log) << "Holding back completed task " << *holdback;
        
        writable_tasks.pop_back();
        
        write(std::move(writable_tasks), temp_vcf);
    } else {
        if (debug_log) stream(*debug_log) << "Buffering completed task " << task;
        buffered_tasks.emplace(contig_region(task), std::move(task));
    }
}

using FutureCompletedTasks = std::vector<std::future<CompletedTask>>;

using RemainingTaskMap = std::map<ContigName, std::deque<CompletedTask>>;

auto extract_remaining_tasks(FutureCompletedTasks& futures, CompletedTaskMap& buffered_tasks)
{
    std::deque<CompletedTask> tasks {};
    
    std::transform(std::begin(futures), std::end(futures), std::back_inserter(tasks),
                   [] (auto& fut) {
                       try {
                           return fut.get();
                       } catch (...) {
                           throw; // TODO: which exceptions can we recover from?
                       }
                   });
    
    futures.clear();
    futures.shrink_to_fit();
    
    for (auto& p : buffered_tasks) {
        std::transform(std::make_move_iterator(std::begin(p.second)),
                       std::make_move_iterator(std::end(p.second)),
                       std::back_inserter(tasks),
                       [] (auto&& p) { return std::move(p.second); });
        p.second.clear();
    }
    
    buffered_tasks.clear();
    
    std::sort(std::begin(tasks), std::end(tasks),
              [] (const auto& lhs, const auto& rhs) {
                  return contig_region(lhs) < contig_region(rhs);
              });
    
    RemainingTaskMap result {};
    
    for (auto&& task : tasks) {
        result[contig_name(task)].push_back(std::move(task));
    }
    
    return result;
}

void resolve_connecting_calls(RemainingTaskMap& remaining_tasks,
                              const ContigCallingComponentFactoryMap& calling_components)
{
    for (auto& p : remaining_tasks) {
        resolve_connecting_calls(p.second, calling_components.at(p.first));
    }
}

void write(RemainingTaskMap&& remaining_tasks, TempVcfWriterMap& temp_vcfs)
{
    for (auto& p : remaining_tasks) {
        write(std::move(p.second), temp_vcfs.at(p.first));
    }
}

void write_remaining_tasks(FutureCompletedTasks& futures, CompletedTaskMap& buffered_tasks,
                           TempVcfWriterMap& temp_vcfs,
                           const ContigCallingComponentFactoryMap& calling_components)
{
    auto remaining_tasks = extract_remaining_tasks(futures, buffered_tasks);
    resolve_connecting_calls(remaining_tasks, calling_components);
    write(std::move(remaining_tasks), temp_vcfs);
}

template <typename K>
auto extract(std::unordered_map<K, VcfWriter>& vcfs)
{
    std::vector<VcfWriter> result {};
    result.reserve(vcfs.size());
    
    for (auto& p : vcfs) result.push_back(std::move(p.second));
    
    vcfs.clear();
    
    return result;
}

template <typename K>
auto extract_as_readers(std::unordered_map<K, VcfWriter>&& vcfs)
{
    auto tmp = extract(vcfs);
    vcfs.clear();
    return writers_to_readers(tmp);
}

template <typename K>
void merge(std::unordered_map<K, VcfWriter>&& temp_vcf_writers, GenomeCallingComponents& components)
{
    auto temp_readers = extract_as_readers(std::move(temp_vcf_writers));
    merge(temp_readers, components.output(), components.contigs_in_output_order());
}

void run_octopus_multi_threaded(GenomeCallingComponents& components)
{
    static auto debug_log = get_debug_log();
    
    const auto num_task_threads = calculate_num_task_threads(components);
    
    // TODO: start running tasks as soon as they are added into the task map,
    // this will require some additional syncronisation on the task map
    
    auto pending_tasks = make_tasks(components, num_task_threads);
    
    auto temp_vcfs = make_temp_vcf_writers(components);
    
    FutureCompletedTasks futures(num_task_threads);
    
    TaskMap running_tasks {ContigOrder {components.contigs_in_output_order()}};
    CompletedTaskMap buffered_tasks {};
    std::map<ContigName, HoldbackTask> holdbacks {};
    
    // populate the all the maps first so we can make unchecked acceses
    for (const auto& contig : components.contigs_in_output_order()) {
        running_tasks.insert({contig, {}});
        buffered_tasks.insert({contig, {}});
        holdbacks.emplace(contig, boost::none);
    }
    
    SyncPacket sync {};
    
    const auto calling_components = make_contig_calling_component_factory_map(components);
    
    components.progress_meter().start();
    
    while (!pending_tasks.empty()) {
        for (auto& future : futures) {
            if (is_ready(future)) {
                std::lock_guard<std::mutex> lk {sync.mutex};
                
                try {
                    auto completed_task = future.get();
                    
                    const auto& contig = contig_name(completed_task);
                    
                    write_or_buffer(std::move(completed_task), buffered_tasks.at(contig),
                                    running_tasks.at(contig), holdbacks.at(contig),
                                    temp_vcfs.at(contig), calling_components.at(contig));
                } catch (...) {
                    throw; // TODO: which exceptions can we recover from?
                }
                
                --sync.count_finished;
            }
            
            if (!pending_tasks.empty() && !future.valid()) {
                std::lock_guard<std::mutex> lk {sync.mutex};
                
                auto task = pop(pending_tasks);
                
                future = run(task, calling_components.at(contig_name(task))(), sync);
                
                running_tasks.at(contig_name(task)).push(std::move(task));
            }
        }
        
        if (sync.count_finished == 0) {
            std::unique_lock<std::mutex> lk {sync.mutex};
            while (sync.count_finished == 0) { // for spurious wakeup
                sync.cv.wait(lk);
            }
        }
    }
    
    running_tasks.clear();
    holdbacks.clear();
    
    write_remaining_tasks(futures, buffered_tasks, temp_vcfs, calling_components);
    
    components.progress_meter().stop();
    
    merge(std::move(temp_vcfs), components);
}
} // namespace

bool is_multithreaded(const GenomeCallingComponents& components)
{
    return !components.num_threads() || *components.num_threads() > 1;
}

auto make_filter_read_pipe(const GenomeCallingComponents& components)
{
    using std::make_unique;
    
    using namespace readpipe;
    
    ReadTransformer transformer {};
    transformer.register_transform(MaskSoftClipped {});
    
    using ReadFilterer = ReadPipe::ReadFilterer;
    
    ReadFilterer filterer {};
    filterer.add(make_unique<HasValidQualities>());
    filterer.add(make_unique<HasWellFormedCigar>());
    filterer.add(make_unique<IsMapped>());
    filterer.add(make_unique<IsNotMarkedQcFail>());
    filterer.add(make_unique<IsNotMarkedDuplicate>());
    filterer.add(make_unique<IsNotDuplicate<ReadFilterer::ReadIterator>>());
    filterer.add(make_unique<IsProperTemplate>());
    
    return ReadPipe {
        components.read_manager(), std::move(transformer), std::move(filterer),
        boost::none, components.samples()
    };
}

auto get_identified_path(const fs::path& base, const std::string& identifier)
{
    const auto old_stem  = base.stem();
    const auto extension = base.extension();
    
    fs::path new_stem;
    
    if (extension.string() == ".gz") {
        new_stem = old_stem.stem().string() + "." + identifier
        + old_stem.extension().string() + extension.string();
    } else {
        new_stem = old_stem.string() + "." + identifier + extension.string();
    }
    
    return base.parent_path() / new_stem;
}

auto get_filtered_path(const fs::path& unfiltered_output_path)
{
    return get_identified_path(unfiltered_output_path, "filtered");
}

auto get_filtered_path(const GenomeCallingComponents& components)
{
    return get_filtered_path(components.output().path());
}

void filter_calls(const GenomeCallingComponents& components)
{
//    assert(!components.output().is_open());
//    
//    const VcfReader calls {components.output().path()};
//    
//    const auto filtered_path = get_filtered_path(components);
//    
//    VcfWriter filtered_calls {filtered_path};
//    
//    const auto read_pipe = make_filter_read_pipe(components);
//    
//    using csr::ThresholdVariantCallFilter;
//    using csr::Measure;
//    using csr::MeasureWrapper;
//    
//    std::vector<MeasureWrapper> measures {};
//    
//    measures.push_back(csr::make_wrapped_measure<csr::QualityByDepth>());
//    
//    auto call_filter = std::make_unique<ThresholdVariantCallFilter>(components.reference(),
//                                                                    read_pipe,
//                                                                    std::move(measures),
//                                                                    components.read_buffer_size());
//    
//    call_filter->filter(calls, filtered_calls);
}

auto get_legacy_path(const fs::path& native)
{
    return get_identified_path(native, "legacy");
}
        
            void print_foo(std::ostream& os)
            {
                os << "foo" << "\n";
            }

void run_octopus(OptionMap& options)
{
    DEBUG_MODE = options::is_debug_mode(options);
    TRACE_MODE = options::is_trace_mode(options);
    
    static auto debug_log = get_debug_log();
    
    logging::InfoLogger info_log {};
    
    const auto start = std::chrono::system_clock::now();
    
    auto end = start;
    
    std::size_t search_size {0};
    
    using utils::TimeInterval;
    
    // open scope to ensure calling components are destroyed before end message
    {
        auto components = collate_genome_calling_components(options);
        
        if (!components) return;
        
        end = std::chrono::system_clock::now();
        
        stream(info_log) << "Done initialising calling components in " << TimeInterval {start, end};
        
        log_startup_info(*components);
        
        if (debug_log) {
            print_input_regions(stream(*debug_log), components->search_regions());
        }
        
        write_caller_output_header(*components, options::call_sites_only(options));
        
        //options.clear();
        
        try {
            if (is_multithreaded(*components)) {
                if (DEBUG_MODE) {
                    logging::WarningLogger warn_log {};
                    warn_log << "Running in parallel mode can make debug log difficult to interpret";
                }
                run_octopus_multi_threaded(*components);
            } else {
                run_octopus_single_threaded(*components);
            }
            
            components->output().close();
        } catch (const std::exception& e) {
            try {
                if (debug_log) {
                    *debug_log << "Encounted an error, attempting to cleanuo";
                }
                cleanup(*components);
            } catch (...) {}
            throw;
        }
        
        if (!options.at("disable-call-filtering").as<bool>()) {
            filter_calls(*components);
            
            if (options.at("legacy").as<bool>()) {
                const auto filtered_path = get_filtered_path(*components);
                
                const VcfReader native {filtered_path};
                VcfWriter legacy {get_legacy_path(native.path())};
                
                convert_to_legacy(native, legacy);
            }
        } else if (options.at("legacy").as<bool>()) {
            const VcfReader native {components->output().path()};
            VcfWriter legacy {get_legacy_path(native.path())};
            
            convert_to_legacy(native, legacy);
        }
        
        cleanup(*components);
        
        end = std::chrono::system_clock::now();
        
        search_size = sum_region_sizes(components->search_regions());
    }
    
    stream(info_log) << "Finished processing " << search_size << "bp, total runtime " << TimeInterval {start, end};
}
} // namespace octopus


