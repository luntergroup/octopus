// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "octopus.hpp"

#include <vector>
#include <deque>
#include <queue>
#include <map>
#include <set>
#include <algorithm>
#include <numeric>
#include <memory>
#include <functional>
#include <cstddef>
#include <typeinfo>
#include <thread>
#include <future>
#include <condition_variable>
#include <mutex>
#include <atomic>
#include <chrono>
#include <sstream>
#include <iostream>
#include <cassert>

#include <boost/optional.hpp>

#include "date/tz.h"
#include "date/ptz.h"

#include "config/common.hpp"
#include "basics/genomic_region.hpp"
#include "basics/ploidy_map.hpp"
#include "concepts/mappable.hpp"
#include "containers/mappable_flat_multi_set.hpp"
#include "containers/mappable_map.hpp"
#include "io/reference/reference_genome.hpp"
#include "io/read/read_manager.hpp"
#include "readpipe/read_pipe_fwd.hpp"
#include "readpipe/buffered_read_pipe.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/read_stats.hpp"
#include "utils/append.hpp"
#include "config/octopus_vcf.hpp"
#include "core/callers/caller_factory.hpp"
#include "core/callers/caller.hpp"
#include "utils/maths.hpp"
#include "logging/progress_meter.hpp"
#include "logging/logging.hpp"
#include "logging/error_handler.hpp"
#include "core/tools/vcf_header_factory.hpp"
#include "io/variant/vcf.hpp"
#include "utils/timing.hpp"
#include "exceptions/program_error.hpp"
#include "exceptions/system_error.hpp"
#include "csr/filters/variant_call_filter.hpp"
#include "csr/filters/variant_call_filter_factory.hpp"
#include "readpipe/buffered_read_pipe.hpp"
#include "core/tools/bam_realigner.hpp"
#include "core/tools/indel_profiler.hpp"
#include "utils/thread_pool.hpp"

#include "timers.hpp" // BENCHMARK

namespace octopus {

using logging::get_debug_log;

namespace {

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

bool apply_csr(const GenomeCallingComponents& components) noexcept
{
    return static_cast<bool>(components.filtered_output());
}

using CallTypeSet = std::set<std::type_index>;

template <class Duration>
inline
auto
make_zoned(const Posix::time_zone& tz, const date::sys_time<Duration>& st)
{
    return date::zoned_time<typename std::common_type<Duration, std::chrono::seconds>::type,
                            Posix::time_zone>{tz, st};
}

std::string get_current_time_str()
{
    using namespace date;
    using namespace std::chrono;
    std::ostringstream ss {};
    try {
        const auto time = make_zoned(current_zone(), system_clock::now());
        ss << time;
    } catch (const std::exception& e) {
        // https://github.com/HowardHinnant/date/issues/310
        static bool warned {false};
        if (!warned) {
            logging::WarningLogger warn_log {};
            warn_log << "Failed to find system time zone information. Assuming UTC.";
            warned = true;
        }
        const auto time = make_zoned(Posix::time_zone{"UTC0"}, system_clock::now());
        ss << time;
    }
    return ss.str();
}

VcfHeader make_vcf_header(const std::vector<SampleName>& samples,
                          const std::vector<GenomicRegion::ContigName>& contigs,
                          const ReferenceGenome& reference,
                          const CallTypeSet& call_types,
                          const UserCommandInfo& info)
{
    auto builder = vcf::make_header_template().set_samples(samples);
    for (const auto& contig : contigs) {
        builder.add_contig(contig, {{"length", std::to_string(reference.contig_size(contig))}});
    }
    builder.add_basic_field("reference", reference.name());
    builder.add_structured_field("octopus", {
        {"version", to_string(config::Version, false)},
        {"command", '"' + info.command + '"'},
        {"options", '"' + info.options + '"'},
        {"date", '"' + get_current_time_str() + '"'}
    });
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
                          const CallTypeSet& call_types,
                          const UserCommandInfo& info)
{
    return make_vcf_header(samples, std::vector<GenomicRegion::ContigName> {contig},
                           reference, call_types, info);
}

VcfHeader read_vcf_header(const VcfReader::Path& vcf_filename)
{
    VcfReader vcf {vcf_filename};
    return vcf.fetch_header();
}

VcfHeader make_vcf_header(const GenomeCallingComponents::Path& filter_vcf,
                          const UserCommandInfo& info)
{
    VcfHeader::Builder builder {read_vcf_header(filter_vcf)};
    builder.add_structured_field("octopus", {
        {"version", to_string(config::Version, false)},
        {"command", '"' + info.command + '"'},
        {"options", '"' + info.options + '"'},
        {"date", '"' + get_current_time_str() + '"'}
    });
    return builder.build_once();
}

bool has_reads(const GenomicRegion& region, ContigCallingComponents& components)
{
    return components.read_manager.get().has_reads(components.samples.get(), region);
}

auto get_call_types(const GenomeCallingComponents& components, const std::vector<ContigName>& contigs)
{
    CallTypeSet result {};
    for (const auto& contig : contigs) {
        const auto tmp_caller = components.caller_factory().make(contig);
        auto caller_call_types = tmp_caller->call_types();
        result.insert(std::begin(caller_call_types), std::end(caller_call_types));
    }
    return result;
}

void write_caller_output_header(GenomeCallingComponents& components, const UserCommandInfo& info)
{
    const auto call_types = get_call_types(components, components.contigs());
    if (components.sites_only() && !apply_csr(components)) {
        components.output() << make_vcf_header({}, components.contigs(), components.reference(),
                                               call_types, info);
    } else if (components.filter_request()) {
        components.output() << make_vcf_header(*components.filter_request(), info);
    } else {
        components.output() << make_vcf_header(components.samples(), components.contigs(),
                                               components.reference(), call_types, info);
    }
}

std::string get_caller_name(const GenomeCallingComponents& components)
{
    assert(!components.contigs().empty());
    const auto& test_contig = components.contigs().front();
    const auto test_caller = components.caller_factory().make(test_contig);
    return test_caller->name(); // There can only be one caller type per run
}

auto my_hardware_concurrency()
{
    // From https://stackoverflow.com/a/31362324/2970186
    std::ifstream cpuinfo {"/proc/cpuinfo"};
    return std::count(std::istream_iterator<std::string> {cpuinfo},
                      std::istream_iterator<std::string> {},
                      std::string {"processor"});
}

unsigned hardware_concurrency()
{
    unsigned int cores = std::thread::hardware_concurrency();
    return cores ? cores : my_hardware_concurrency();
}

class UnknownHardwareConcurrency : public SystemError
{
    std::string do_where() const override { return "hardware_concurrency"; }
    std::string do_why() const override { return "Unable to detect the number of hardware threads"; };
    std::string do_help() const override { return "Explicitly set the number of threads in the command line input"; };
};

void log_startup_info(const GenomeCallingComponents& components)
{
    logging::InfoLogger log {};
    std::ostringstream ss {};
    if (!components.samples().empty()) {
        const auto num_samples = components.samples().size();
        if (num_samples == 1) {
            ss << "Detected 1 sample: ";
        } else {
            ss << "Detected " << num_samples << " samples: ";
        }
    }
    std::transform(std::cbegin(components.samples()), std::cend(components.samples()),
                   std::ostream_iterator<std::string> {ss, " "},
                   [] (const auto& sample) -> std::string {
                       return "\"" + sample + "\"";
                   });
    auto str = ss.str();
    str.pop_back(); // the extra whitespace
    log << str;
    stream(log) << "Invoked calling model: " << get_caller_name(components);
    {
        const auto search_size = utils::format_with_commas(sum_region_sizes(components.search_regions()));
        const auto num_threads = components.num_threads();
        auto ls = stream(log);
        if (num_threads) {
            if (*num_threads == 1) {
                ls << "Processing " << search_size << "bp with a single thread";
            } else {
                ls << "Processing " << search_size << "bp with " << *num_threads << " threads";
            }
        } else {
            ls << "Processing " << search_size << "bp with automatic thread management";
        }
        const auto cores = hardware_concurrency();
        if (cores > 0) {
            ls << " (" << cores << " cores detected)";
        }
    }
    auto sl = stream(log);
    auto output_path = components.output().path();
    if (apply_csr(components)) {
        sl << "Writing filtered calls to ";
        output_path = components.filtered_output()->path();
    } else {
        sl << "Writing unfiltered calls to ";
    }
    if (output_path) {
        sl << *output_path;
    } else {
        sl << "stdout";
    }
}

VcfWriter& get_final_output(GenomeCallingComponents& components)
{
    if (apply_csr(components)) {
        return *components.filtered_output();
    } else {
        return components.output();
    }
}

auto get_final_output_path(const GenomeCallingComponents& components)
{
    if (components.filtered_output()) {
        return components.filtered_output()->path();
    } else {
        return components.output().path();
    }
}

void log_finish_info(const GenomeCallingComponents& components, const utils::TimeInterval run_duration)
{
    using utils::TimeInterval;
    logging::InfoLogger info_log {};
    const auto search_size = sum_region_sizes(components.search_regions());
    stream(info_log) << "Finished calling "
                     << utils::format_with_commas(search_size) << "bp, total runtime "
                     << run_duration;
    const auto output_path = get_final_output_path(components);
    if (output_path) stream(info_log) << "Calls have been written to " << *output_path;
}

void write_calls(std::deque<VcfRecord>&& calls, VcfWriter& out)
{
    if (calls.empty()) return;
    static auto debug_log = get_debug_log();
    if (debug_log) stream(*debug_log) << "Writing " << calls.size() << " calls to output";
    const bool was_closed {!out.is_open()};
    if (was_closed) out.open();
    write(calls, out);
    if (was_closed) out.close();
    calls.clear();
    calls.shrink_to_fit();
}

struct WindowConfig
{
    boost::optional<GenomicRegion::Size> min_size = boost::none, max_size = boost::none;
};

static const WindowConfig default_window_config {5'000, 25'000'000};

auto find_max_window(const ContigCallingComponents& components,
                     const GenomicRegion& target_region)
{
    const auto& rm = components.read_manager.get();
    if (!rm.has_reads(components.samples.get(), target_region)) {
        return target_region;
    }
    auto result = rm.find_covered_subregion(components.samples, target_region, components.read_buffer_size);
    if (ends_before(result, target_region)) {
        auto rest = right_overhang_region(target_region, result);
        if (!rm.has_reads(components.samples.get(), rest)) {
            result = target_region;
        }
    }
    return result;
}

auto propose_call_subregion(const ContigCallingComponents& components,
                            const GenomicRegion& remaining_call_region,
                            const WindowConfig& config)
{
    if (is_empty(remaining_call_region)) {
        return remaining_call_region;
    }
    if (config.min_size && size(remaining_call_region) <= *config.min_size) {
        return remaining_call_region;
    }
    auto target = remaining_call_region;
    if (config.max_size && size(target) > *config.max_size) {
        target = head_region(target, *config.max_size);
    }
    const auto max_window = find_max_window(components, target);
    if (ends_before(remaining_call_region, max_window)) {
        return remaining_call_region;
    }
    if (config.min_size && size(max_window) < *config.min_size) {
        return expand_rhs(head_region(max_window), *config.min_size);
    }
    return max_window;
}

auto propose_call_subregion(const ContigCallingComponents& components,
                            const GenomicRegion& current_subregion,
                            const GenomicRegion& input_region,
                            const WindowConfig& config)
{
    assert(contains(input_region, current_subregion));
    return propose_call_subregion(components, right_overhang_region(input_region, current_subregion), config);
}

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
        std::vector<VcfRecord> new_connecting_calls {};
        const auto old_connecting_calls_region = encompassing_region(old_connecting_calls);
        buffer_connecting_calls(old_connecting_calls_region, calls, new_connecting_calls);
        std::deque<VcfRecord> merged_calls {};
        std::set_union(begin(old_connecting_calls), end(old_connecting_calls),
                       begin(new_connecting_calls), end(new_connecting_calls),
                       std::back_inserter(merged_calls));
        old_connecting_calls.clear();
        old_connecting_calls.shrink_to_fit();
        new_connecting_calls.clear();
        new_connecting_calls.shrink_to_fit();
        if (is_consistent(merged_calls)) {
            calls.insert(begin(calls),
                         make_move_iterator(begin(merged_calls)),
                         make_move_iterator(end(merged_calls)));
        } else {
            const auto unresolved_region = encompassing_region(merged_calls);
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
    
    const auto window_config = default_window_config;
    
    std::deque<VcfRecord> calls;
    std::vector<VcfRecord> connecting_calls {};
    auto input_region = components.regions.front();
    auto subregion    = propose_call_subregion(components, input_region, window_config);
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
        
        auto next_subregion = propose_call_subregion(components, subregion, input_region, window_config);
        
        if (is_empty(next_subregion)) {
            ++first_input_region;
            if (first_input_region != last_input_region) {
                input_region = *first_input_region;
                next_subregion = propose_call_subregion(components, input_region, window_config);
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
}

void run_octopus_single_threaded(GenomeCallingComponents& components)
{
    #ifdef BENCHMARK
    init_timers();
    #endif
    components.progress_meter().start();
    for (const auto& contig : components.contigs()) {
        run_octopus_on_contig(ContigCallingComponents {contig, components});
    }
    components.progress_meter().stop();
    #ifdef BENCHMARK
        print_all_timers();
    #endif
}

bool can_use_temp_bcf(const GenomicRegion& region)
{
    // htslib cannot parse ':' in contig names.
    // See https://github.com/samtools/htslib/pull/708
    return std::find(std::cbegin(region.contig_name()), std::cend(region.contig_name()), ':') == std::cend(region.contig_name());
}

auto create_unique_temp_output_file_path(const GenomicRegion& region,
                                         const GenomeCallingComponents& components)
{
    auto result = *components.temp_directory();
    const auto begin   = std::to_string(region.begin());
    const auto end     = std::to_string(region.end());
    boost::filesystem::path file_name {region.contig_name() + "_" + begin + "-" + end + "_temp"};
    
    // Hack for htslib ':' parsing issues
    if (can_use_temp_bcf(region)) {
        file_name += ".bcf";
    } else {
        file_name += ".vcf";
    }
    
    result /= file_name;
    return result;
}

VcfHeader make_temp_vcf_header(const GenomeCallingComponents& components, const GenomicRegion& region)
{
    const auto call_types = get_call_types(components, {region.contig_name()});
    return make_vcf_header(components.samples(), region.contig_name(), components.reference(), call_types, {"octopus-internal", ""});
}

VcfWriter create_unique_temp_output_file(const GenomicRegion& region, const GenomeCallingComponents& components)
{
    return {create_unique_temp_output_file_path(region, components), make_temp_vcf_header(components, region)};
}

VcfWriter create_unique_temp_output_file(const GenomicRegion::ContigName& contig, const GenomeCallingComponents& components)
{
    return create_unique_temp_output_file(components.reference().contig_region(contig), components);
}

using TempVcfWriterMap = std::unordered_map<ContigName, VcfWriter>;

TempVcfWriterMap make_temp_vcf_writers(const GenomeCallingComponents& components)
{
    if (!components.temp_directory()) {
        throw std::runtime_error {"Could not make temp writers"};
    }
    TempVcfWriterMap result {};
    result.reserve(components.contigs().size());
    for (const auto& contig : components.contigs()) {
        auto contig_writer = create_unique_temp_output_file(contig, components);
        contig_writer.close();
        result.emplace(contig, std::move(contig_writer));
    }
    return result;
}

struct Task : public Mappable<Task>
{
    GenomicRegion region;
    ExecutionPolicy policy;
    
    Task() = delete;
    
    Task(GenomicRegion region, ExecutionPolicy policy = ExecutionPolicy::seq)
    : region {std::move(region)}
    , policy {policy}
    {};
    
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
    : contigs_ {std::cbegin(contigs)
    , std::cend(contigs)}
    {}
    
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

auto count_tasks(const TaskMap& tasks) noexcept
{
    return std::accumulate(std::cbegin(tasks), std::cend(tasks), std::size_t {0},
                           [] (auto curr, const auto& p) noexcept { return curr + p.second.size(); });
}

struct TaskMakerSyncPacket
{
    TaskMakerSyncPacket() : batch_size_hint {1}, waiting {true}, num_tasks {0}, finished {}, all_done {false} {}
    std::condition_variable cv;
    std::mutex mutex;
    bool ready = true;
    std::atomic_uint batch_size_hint; // Only read by task maker
    std::atomic_bool waiting; // Only read by task maker
    std::atomic_uint num_tasks;
    std::unordered_map<ContigName, bool> finished;
    std::atomic_bool all_done;
};

void make_region_tasks(const GenomicRegion& region,
                       const ContigCallingComponents& components,
                       const ExecutionPolicy policy,
                       TaskQueue& result,
                       TaskMakerSyncPacket& sync,
                       const bool last_region_in_contig,
                       const bool last_contig,
                       const WindowConfig& window_config)
{
    std::unique_lock<std::mutex> lock {sync.mutex, std::defer_lock};
    auto subregion = propose_call_subregion(components, region, window_config);
    if (ends_equal(subregion, region)) {
        lock.lock();
        sync.cv.wait(lock, [&] () { return sync.ready; });
        result.emplace(std::move(subregion), policy);
        ++sync.num_tasks;
        if (last_region_in_contig) {
            sync.finished.at(region.contig_name()) = true;
            if (last_contig) sync.all_done = true;
        }
        lock.unlock();
        sync.cv.notify_one();
    } else {
        std::deque<GenomicRegion> batch {};
        batch.push_back(subregion);
        bool done {false};
        while (true) {
            while (batch.size() < std::max(sync.batch_size_hint.load(), 1u) || !sync.waiting) {
                subregion = propose_call_subregion(components, subregion, region, window_config);
                batch.push_back(subregion);
                assert(!ends_before(region, subregion));
                if (ends_equal(subregion, region)) {
                    done = true;
                    break;
                }
            }
            assert(!batch.empty());
            assert(!lock.owns_lock());
            lock.lock();
            sync.cv.wait(lock, [&] () { return sync.ready; });
            for (auto&& r : batch) result.emplace(std::move(r), policy);
            sync.num_tasks += batch.size();
            if (done) {
                if (last_region_in_contig) {
                    sync.finished.at(region.contig_name()) = true;
                    if (last_contig) sync.all_done = true;
                } else {
                    assert(!last_contig);
                }
                lock.unlock();
                sync.cv.notify_one();
                break;
            } else {
                batch.clear();
                lock.unlock();
                sync.cv.notify_one();
            }
        }
    }
}

void make_contig_tasks(const ContigCallingComponents& components,
                       const ExecutionPolicy policy,
                       TaskQueue& result,
                       TaskMakerSyncPacket& sync,
                       const bool last_contig,
                       const WindowConfig& window_config)
{
    if (components.regions.empty()) return;
    std::for_each(std::cbegin(components.regions), std::prev(std::cend(components.regions)), [&] (const auto& region) {
        make_region_tasks(region, components, policy, result, sync, false, last_contig, window_config);
    });
    make_region_tasks(components.regions.back(), components, policy, result, sync, true, last_contig, window_config);
}

ExecutionPolicy make_execution_policy(const GenomeCallingComponents& components)
{
    if (components.num_threads()) {
        return ExecutionPolicy::seq;
    }
    return ExecutionPolicy::par;
}

auto make_contig_components(const ContigName& contig, GenomeCallingComponents& components, const unsigned num_threads)
{
    ContigCallingComponents result {contig, components};
    result.read_buffer_size /= num_threads;
    return result;
}

void make_tasks_helper(TaskMap& tasks,
                       std::vector<ContigName> contigs,
                       GenomeCallingComponents& components,
                       const unsigned num_threads,
                       ExecutionPolicy execution_policy,
                       TaskMakerSyncPacket& sync)
{
    const auto window_config = default_window_config;
    try {
        static auto debug_log = get_debug_log();
        if (debug_log) stream(*debug_log) << "Making tasks for " << contigs.size() << " contigs";
        for (std::size_t i {0}; i < contigs.size(); ++i) {
            const auto& contig = contigs[i];
            if (debug_log) stream(*debug_log) << "Making tasks for contig " << contig;
            auto contig_components = make_contig_components(contig, components, num_threads);
            make_contig_tasks(contig_components, execution_policy, tasks[contig], sync, i == contigs.size() - 1, window_config);
            if (debug_log) stream(*debug_log) << "Finished making tasks for contig " << contig;
        }
        if (debug_log) *debug_log << "Finished making tasks";
    } catch (const Error& e) {
        log_error(e);
        logging::FatalLogger fatal_log {};
        fatal_log << "Encountered error in task maker thread. Calling terminate";
        std::terminate();
    } catch (const std::exception& e) {
        log_error(e);
        logging::FatalLogger fatal_log {};
        fatal_log << "Encountered error in task maker thread. Calling terminate";
        std::terminate();
    } catch (...) {
        logging::FatalLogger fatal_log {};
        fatal_log << "Encountered error in task maker thread. Calling terminate";
        std::terminate();
    }
}

std::thread
make_task_maker_thread(TaskMap& tasks,
                       GenomeCallingComponents& components,
                       const unsigned num_threads,
                       TaskMakerSyncPacket& sync)
{
    auto contigs = components.contigs();
    if (contigs.empty()) {
        sync.all_done = true;
        return std::thread {};
    }
    sync.finished.reserve(contigs.size());
    for (const auto& contig : contigs) {
        sync.finished.emplace(contig, false);
    }
    return std::thread {make_tasks_helper, std::ref(tasks), std::move(contigs), std::ref(components),
                        num_threads, make_execution_policy(components), std::ref(sync)};
}

unsigned calculate_num_task_threads(const GenomeCallingComponents& components)
{
    if (components.num_threads()) {
        return *components.num_threads();
    }
    const auto num_cores = hardware_concurrency();
    if (num_cores == 0) {
        throw UnknownHardwareConcurrency {};
    }
    return num_cores;
}

Task pop(TaskMap& tasks, TaskMakerSyncPacket& sync)
{
    assert(!tasks.empty());
    std::unique_lock<std::mutex> lock {sync.mutex};
    sync.ready = false;
    assert(sync.num_tasks > 0);
    const auto contig_task_itr = std::begin(tasks);
    assert(!contig_task_itr->second.empty());
    const auto result = std::move(contig_task_itr->second.front());
    contig_task_itr->second.pop();
    if (sync.finished.at(contig_task_itr->first) && contig_task_itr->second.empty()) {
        static auto debug_log = get_debug_log();
        if (debug_log) stream(*debug_log) << "Finished calling contig " << contig_task_itr->first;
        tasks.erase(contig_task_itr);
    }
    --sync.num_tasks;
    sync.ready = true;
    lock.unlock();
    sync.cv.notify_one();
    return result;
}

struct CompletedTask : public Task
{
    CompletedTask(Task task) : Task {std::move(task)}, calls {}, runtime {} {}
    std::deque<VcfRecord> calls;
    utils::TimeInterval runtime;
};

std::string duration(const CompletedTask& task)
{
    std::ostringstream ss {};
    ss << task.runtime;
    return ss.str();
}

std::ostream& operator<<(std::ostream& os, const CompletedTask& task)
{
    os << task.region;
    return os;
}

struct CallerSyncPacket
{
    CallerSyncPacket() : num_finished {0} {}
    std::condition_variable cv;
    std::mutex mutex;
    std::atomic_uint num_finished;
};

template<typename R>
bool is_ready(const std::future<R>& f)
{
    return f.valid() && f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
}

auto run(Task task, ContigCallingComponents components, CallerSyncPacket& sync, ThreadPool& workers)
{
    static auto debug_log = get_debug_log();
    if (debug_log) stream(*debug_log) << "Spawning task " << task;
    return workers.push([task = std::move(task), components = std::move(components), &sync, &workers] () {
        try {
            CompletedTask result {task};
            result.runtime.start = std::chrono::system_clock::now();
            result.calls = components.caller->call(task.region, components.progress_meter, workers);
            result.runtime.end = std::chrono::system_clock::now();
            std::unique_lock<std::mutex> lock {sync.mutex};
            ++sync.num_finished;
            lock.unlock();
            sync.cv.notify_all();
            return result;
        } catch (const std::exception& e) {
            logging::ErrorLogger error_log {};
            stream(error_log) << "Encountered a problem whilst calling " << task << "(" << e.what() << ")";
            using namespace std::chrono_literals;
            std::this_thread::sleep_for(2s); // Try to make sure the error is logged before raising
            std::unique_lock<std::mutex> lock {sync.mutex};
            ++sync.num_finished;
            lock.unlock();
            sync.cv.notify_all();
            throw;
        }
    });
}

using CompletedTaskMap = std::map<ContigName, std::map<ContigRegion, CompletedTask>>;
using HoldbackTask = boost::optional<std::reference_wrapper<const CompletedTask>>;

auto get_writable_completed_tasks(CompletedTask&& task, CompletedTaskMap::mapped_type& buffered_tasks,
                                  TaskQueue& running_tasks, HoldbackTask& holdback)
{
    std::deque<CompletedTask> result {std::move(task)};
    while (!running_tasks.empty()) {
        const auto itr = buffered_tasks.find(contig_region(running_tasks.front()));
        if (itr != std::end(buffered_tasks)) {
            result.push_back(std::move(itr->second));
            buffered_tasks.erase(itr);
            running_tasks.pop();
        } else {
            break;
        }
    }
    if (holdback) {
        const auto itr = buffered_tasks.find(contig_region(holdback->get()));
        assert(itr->second < result.front());
        result.push_front(std::move(itr->second));
        buffered_tasks.erase(itr);
        holdback = boost::none;
    }
    return result;
}

using ContigCallingComponentFactory    = std::function<ContigCallingComponents()>;
using ContigCallingComponentFactoryMap = std::map<ContigName, ContigCallingComponentFactory>;

auto make_contig_calling_component_factory_map(GenomeCallingComponents& components)
{
    ContigCallingComponentFactoryMap result {};
    for (const auto& contig : components.contigs()) {
        result.emplace(contig, [&components, contig] () -> ContigCallingComponents
                       { return ContigCallingComponents {contig, components}; });
    }
    return result;
}

auto find_first_lhs_connecting(const std::deque<VcfRecord>& lhs_calls, const GenomicRegion& rhs_region)
{
    const auto rhs_begin = mapped_begin(rhs_region);
    return std::find_if(std::cbegin(lhs_calls), std::cend(lhs_calls),
                        [&rhs_begin] (const auto& call) { return mapped_end(call) > rhs_begin; });
}

auto find_last_rhs_connecting(const GenomicRegion& lhs_region, const std::deque<VcfRecord>& rhs_calls)
{
    const auto lhs_end = mapped_end(lhs_region);
    return std::find_if_not(std::cbegin(rhs_calls), std::cend(rhs_calls),
                            [&lhs_end] (const auto& call) { return mapped_begin(call) < lhs_end; });
}

void resolve_connecting_calls(CompletedTask& lhs, CompletedTask& rhs,
                              const ContigCallingComponentFactory& calling_components)
{
    static auto debug_log = get_debug_log();
    if (lhs.calls.empty() || rhs.calls.empty()) return;
    const auto lhs_region = encompassing_region(lhs.calls);
    const auto rhs_region = encompassing_region(rhs.calls);
    const auto first_lhs_connecting = find_first_lhs_connecting(lhs.calls, rhs_region);
    const auto last_rhs_connecting  = find_last_rhs_connecting(lhs_region, rhs.calls);
    if (debug_log) {
        const auto num_lhs_conflicting = std::distance(first_lhs_connecting, std::cend(lhs.calls));
        const auto num_rhs_conflicting = std::distance(std::cbegin(rhs.calls), last_rhs_connecting);
        if (num_lhs_conflicting + num_rhs_conflicting > 0) {
            stream(*debug_log) << "Resolving connecting calls between tasks " 
                           << lhs << "(" << num_lhs_conflicting << " conflicting)"
                           << " & " << rhs << "(" << num_rhs_conflicting << " conflicting)";
        }
    }
    // The general stratergy is to keep RHS variant calls otherwise we might mess up phase sets of downstream calls.
    // However, we don't need to keep RHS reference call before the first variant call, so
    // we should prefer to keep LHS variant calls in this region.
    const static auto is_refcall = [] (const VcfRecord& call) { return call.is_refcall(); };
    const auto first_rhs_variant = std::find_if_not(std::cbegin(rhs.calls), last_rhs_connecting, is_refcall);
    auto first_lhs_remove = first_lhs_connecting;
    if (first_rhs_variant != std::cbegin(rhs.calls)) {
        const auto rhs_ref_region = first_rhs_variant != std::cend(rhs.calls) ?
             closed_region(rhs.calls.front(), *first_rhs_variant) : rhs_region;
        const auto rhs_keep_region = right_overhang_region(rhs_region, rhs_ref_region);
        first_lhs_remove = std::find_if(first_lhs_connecting, std::cend(lhs.calls),
                 [&] (const auto& call) { return overlaps(call, rhs_keep_region); });
        if (first_lhs_remove != std::cbegin(lhs.calls)) {
            const auto lhs_keep_region = encompassing_region(std::cbegin(lhs.calls), first_lhs_remove);
            const auto last_rhs_remove = std::find_if_not(std::cbegin(rhs.calls), first_rhs_variant,
                    [&] (const auto& call) { return overlaps(call, lhs_keep_region); });
            if (last_rhs_remove != std::cbegin(rhs.calls)) {
                const auto last_overlapping_rhs_ref = *std::prev(last_rhs_remove);
                rhs.calls.erase(std::cbegin(rhs.calls), last_rhs_remove);
                if (!contains(lhs_keep_region, last_overlapping_rhs_ref)) {
                    // If the last RHS refcall overlapping with a LHS call we'll keep is a partial
                    // overlap then we should still keep the non-overlapping positions.
                    auto squashed_refcall = VcfRecord::Builder(last_overlapping_rhs_ref).set_pos(mapped_end(lhs_keep_region)).build_once();
                    rhs.calls.push_front(std::move(squashed_refcall));
                }
            }
        }
    }
    lhs.calls.erase(first_lhs_remove, std::cend(lhs.calls));
}

void resolve_connecting_calls(std::deque<CompletedTask>& adjacent_tasks,
                              const ContigCallingComponentFactory& calling_components)
{
    if (adjacent_tasks.size() < 2) return;
    assert(std::is_sorted(std::cbegin(adjacent_tasks), std::cend(adjacent_tasks)));
    auto lhs = std::begin(adjacent_tasks);
    std::for_each(std::next(lhs), std::end(adjacent_tasks),
                  [&] (auto& rhs) { resolve_connecting_calls(*lhs++, rhs, calling_components); });
}

struct TaskWriterSyncPacket
{
    std::condition_variable cv;
    std::mutex mutex;
    std::deque<CompletedTask> tasks = {};
    bool done = false;
    bool completed = false;
};

void write(std::deque<CompletedTask>& tasks, TempVcfWriterMap& writers)
{
    static auto debug_log = get_debug_log();
    for (auto&& task : tasks) {
        if (debug_log) {
            stream(*debug_log) << "Writing completed task " << task << " that finished in " << duration(task);
        }
        auto& writer = writers.at(contig_name(task));
        write_calls(std::move(task.calls), writer);
    }
    tasks.clear();
}

void write_temp_vcf_helper(TempVcfWriterMap& writers, TaskWriterSyncPacket& sync)
{
    static auto debug_log = get_debug_log();
    try {
        std::unique_lock<std::mutex> lock {sync.mutex, std::defer_lock};
        std::deque<CompletedTask> buffer {};
        while (!sync.done) {
            lock.lock();
            sync.cv.wait(lock, [&] () { return !sync.tasks.empty() || sync.done; });
            assert(buffer.empty());
            std::swap(sync.tasks, buffer);
            lock.unlock();
            sync.cv.notify_one();
            write(buffer, writers);
        }
        if (debug_log) *debug_log << "Task writer finished";
        lock.lock();
        sync.completed = true;
        sync.cv.notify_all();
    } catch (const Error& e) {
        log_error(e);
        logging::FatalLogger fatal_log {};
        fatal_log << "Encountered error in task writer thread. Calling terminate";
        std::terminate();
    } catch (const std::exception& e) {
        log_error(e);
        logging::FatalLogger fatal_log {};
        fatal_log << "Encountered error in task writer thread. Calling terminate";
        std::terminate();
    } catch (...) {
        logging::FatalLogger fatal_log {};
        fatal_log << "Encountered error in task writer thread. Calling terminate";
        std::terminate();
    }
}

std::thread make_task_writer_thread(TempVcfWriterMap& temp_writers, TaskWriterSyncPacket& writer_sync)
{
    return std::thread {write_temp_vcf_helper, std::ref(temp_writers), std::ref(writer_sync)};
}

void write(std::deque<CompletedTask>&& tasks, VcfWriter& temp_vcf)
{
    static auto debug_log = get_debug_log();
    for (auto&& task : tasks) {
        if (debug_log) stream(*debug_log) << "Writing completed task " << task << " that finished in " << duration(task);
        write_calls(std::move(task.calls), temp_vcf);
    }
}

void write(std::deque<CompletedTask>&& tasks, TaskWriterSyncPacket& sync)
{
    std::unique_lock<std::mutex> lock {sync.mutex};
    utils::append(std::move(tasks), sync.tasks);
    lock.unlock();
    sync.cv.notify_one();
}

// A CompletedTask can only be written if all proceeding tasks have completed (either written or buffered)
void write_or_buffer(CompletedTask&& task, CompletedTaskMap::mapped_type& buffered_tasks,
                     TaskQueue& running_tasks, HoldbackTask& holdback,
                     TaskWriterSyncPacket& sync, const ContigCallingComponentFactory& calling_components)
{
    static auto debug_log = get_debug_log();
    if (is_same_region(task, running_tasks.front())) {
        running_tasks.pop();
        auto writable_tasks = get_writable_completed_tasks(std::move(task), buffered_tasks, running_tasks, holdback);
        assert(holdback == boost::none);
        resolve_connecting_calls(writable_tasks, calling_components);
        // Keep the last task buffered to enable connection resolution when the next task finishes
        const auto p = buffered_tasks.emplace(contig_region(writable_tasks.back()), std::move(writable_tasks.back()));
        assert(p.second);
        holdback = p.first->second;
        if (debug_log) stream(*debug_log) << "Holding back completed task " << *holdback;
        writable_tasks.pop_back();
        write(std::move(writable_tasks), sync);
    } else {
        if (debug_log) stream(*debug_log) << "Buffering completed task " << task;
        buffered_tasks.emplace(contig_region(task), std::move(task));
    }
}

void wait_until_finished(TaskWriterSyncPacket& sync)
{
    std::unique_lock<std::mutex> lock {sync.mutex};
    sync.cv.wait(lock, [&] () { return sync.tasks.empty(); });
    sync.done = true;
    lock.unlock();
    sync.cv.notify_one();
    lock.lock();
    sync.cv.wait(lock, [&] () { return sync.completed; });
}

using FutureCompletedTasks = std::vector<std::future<CompletedTask>>;
using RemainingTaskMap = std::map<ContigName, std::deque<CompletedTask>>;

void extract_remaining_future_tasks(FutureCompletedTasks& futures, std::deque<CompletedTask>& result)
{
    const auto itr = std::remove_if(std::begin(futures), std::end(futures), [] (const auto& f) { return !f.valid(); });
    std::transform(std::begin(futures), itr, std::back_inserter(result), [] (auto& fut) { return fut.get(); });
    futures.clear();
    futures.shrink_to_fit();
}

void extract_buffered_tasks(CompletedTaskMap& buffered_tasks, std::deque<CompletedTask>& result)
{
    for (auto& p : buffered_tasks) {
        std::transform(std::make_move_iterator(std::begin(p.second)), std::make_move_iterator(std::end(p.second)),
                       std::back_inserter(result), [] (auto&& p) { return std::move(p.second); });
        p.second.clear();
    }
    buffered_tasks.clear();
}

void sort_ignoring_contig_name(std::deque<CompletedTask>& tasks)
{
    std::sort(std::begin(tasks), std::end(tasks),
              [] (const auto& lhs, const auto& rhs) { return contig_region(lhs) < contig_region(rhs); });
}

RemainingTaskMap make_map(std::deque<CompletedTask>& tasks)
{
    sort_ignoring_contig_name(tasks);
    RemainingTaskMap result {};
    for (auto&& task : tasks) {
        result[contig_name(task.region)].push_back(std::move(task));
    }
    return result;
}

RemainingTaskMap extract_remaining_tasks(FutureCompletedTasks& futures, CompletedTaskMap& buffered_tasks)
{
    std::deque<CompletedTask> tasks {};
    extract_remaining_future_tasks(futures, tasks);
    extract_buffered_tasks(buffered_tasks, tasks);
    return make_map(tasks);
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

void write_remaining_tasks(FutureCompletedTasks& futures, CompletedTaskMap& buffered_tasks, TempVcfWriterMap& temp_vcfs,
                           const ContigCallingComponentFactoryMap& calling_components)
{
    static auto debug_log = get_debug_log();
    if (debug_log) stream(*debug_log) << "Waiting for " << futures.size() << " running tasks to finish";
    auto remaining_tasks = extract_remaining_tasks(futures, buffered_tasks);
    resolve_connecting_calls(remaining_tasks, calling_components);
    write(std::move(remaining_tasks), temp_vcfs);
}

auto extract_writers(TempVcfWriterMap&& vcfs)
{
    std::vector<VcfWriter> result {};
    result.reserve(vcfs.size());
    for (auto&& p : vcfs) {
        result.push_back(std::move(p.second));
    }
    vcfs.clear();
    return result;
}

auto extract_as_readers(TempVcfWriterMap&& vcfs)
{
    return writers_to_readers(extract_writers(std::move(vcfs)), false);
}

void merge(TempVcfWriterMap&& temp_vcf_writers, GenomeCallingComponents& components)
{
    static auto debug_log = get_debug_log();
    if (debug_log) stream(*debug_log) << "Merging " << temp_vcf_writers.size() << " temporary VCF files";
    auto temp_readers = extract_as_readers(std::move(temp_vcf_writers));
    merge(temp_readers, components.output(), components.contigs());
}

void run_octopus_multi_threaded(GenomeCallingComponents& components)
{
    using namespace std::chrono_literals;
    static auto debug_log = get_debug_log();
    
    const auto num_task_threads = calculate_num_task_threads(components);
    ThreadPool workers {num_task_threads};
    
    TaskMap pending_tasks {components.contigs()};
    TaskMakerSyncPacket task_maker_sync {};
    task_maker_sync.batch_size_hint = 2 * num_task_threads;
    std::unique_lock<std::mutex> pending_task_lock {task_maker_sync.mutex, std::defer_lock};
    auto task_maker_thread = make_task_maker_thread(pending_tasks, components, num_task_threads, task_maker_sync);
    if (!task_maker_thread.joinable()) {
        logging::FatalLogger fatal_log {};
        fatal_log << "Unable to make task maker thread";
        return;
    }
    task_maker_thread.detach();
    
    FutureCompletedTasks futures(num_task_threads);
    TaskMap running_tasks {ContigOrder {components.contigs()}};
    CompletedTaskMap buffered_tasks {};
    std::map<ContigName, HoldbackTask> holdbacks {};
    // Populate all the maps first so we can make unchecked accesses
    for (const auto& contig : components.contigs()) {
        running_tasks.emplace(contig, TaskMap::mapped_type {});
        buffered_tasks.emplace(contig, CompletedTaskMap::mapped_type {});
        holdbacks.emplace(contig, boost::none);
    }
    
    CallerSyncPacket caller_sync {};
    const auto calling_components = make_contig_calling_component_factory_map(components);
    unsigned num_idle_futures {0};
    
    auto temp_writers = make_temp_vcf_writers(components);
    TaskWriterSyncPacket task_writer_sync {};
    auto task_writer_thread = make_task_writer_thread(temp_writers, task_writer_sync);
    if (!task_writer_thread.joinable()) {
        logging::FatalLogger fatal_log {};
        fatal_log << "Unable to make task writer thread";
        return;
    }
    task_writer_thread.detach();
    
    // Wait for the first task to be made
    const auto tasks_available = [&] () noexcept { return task_maker_sync.num_tasks > 0; };
    while(task_maker_sync.num_tasks == 0) {
        pending_task_lock.lock();
        task_maker_sync.cv.wait(pending_task_lock, tasks_available);
        pending_task_lock.unlock();
    }
    task_maker_sync.batch_size_hint = num_task_threads / 2;
    
    components.progress_meter().start();
    
    while (!task_maker_sync.all_done || task_maker_sync.num_tasks > 0) {
        pending_task_lock.lock();
        assert(count_tasks(pending_tasks) == task_maker_sync.num_tasks);
        if (!task_maker_sync.all_done && task_maker_sync.num_tasks == 0) {
            task_maker_sync.batch_size_hint = std::max(num_idle_futures, num_task_threads / 2);
            if (num_idle_futures < futures.size()) {
                // If there are running futures then it's good periodically check to see if
                // any have finished and process them while we wait for the task maker.
                while (task_maker_sync.num_tasks == 0 && caller_sync.num_finished == 0) {
                    auto now = std::chrono::system_clock::now();
                    task_maker_sync.cv.wait_until(pending_task_lock, now + 5s, tasks_available);
                }
            } else {
                task_maker_sync.cv.wait(pending_task_lock, tasks_available);
            }
        }
        pending_task_lock.unlock();
        num_idle_futures = 0;
        for (auto& future : futures) {
            if (is_ready(future)) {
                auto completed_task = future.get();
                const auto& contig = contig_name(completed_task.region);
                write_or_buffer(std::move(completed_task), buffered_tasks.at(contig),
                                running_tasks.at(contig), holdbacks.at(contig),
                                task_writer_sync, calling_components.at(contig));
                --caller_sync.num_finished;
            }
            if (!future.valid()) {
                pending_task_lock.lock();
                if (task_maker_sync.num_tasks > 0) {
                    pending_task_lock.unlock(); // As pop will need to lock the mutex too == deadlock
                    auto task = pop(pending_tasks, task_maker_sync);
                    future = run(task, calling_components.at(contig_name(task))(), caller_sync, workers);
                    running_tasks.at(contig_name(task)).push(std::move(task));
                } else {
                    pending_task_lock.unlock();
                    ++num_idle_futures;
                }
            }
        }
        // If there are no idle futures then all threads are busy and we must wait for one to finish,
        // otherwise we must have run out of tasks, so we should wait for new ones.
        if (num_idle_futures == 0 && caller_sync.num_finished == 0) {
            task_maker_sync.waiting = false;
            std::unique_lock<std::mutex> lock {caller_sync.mutex};
            caller_sync.cv.wait(lock, [&] () { return caller_sync.num_finished > 0; });
            task_maker_sync.waiting = true;
        } else {
            if (debug_log) stream(*debug_log) << "There are " << num_idle_futures << " idle futures";
        }
    }
    assert(task_maker_sync.num_tasks == 0);
    assert(pending_tasks.empty());
    running_tasks.clear();
    holdbacks.clear(); // holdbacks are just references to buffered tasks
    if (debug_log) *debug_log << "Finished making new tasks. Waiting for task writer to complete existing jobs";
    wait_until_finished(task_writer_sync);
    write_remaining_tasks(futures, buffered_tasks, temp_writers, calling_components);
    components.progress_meter().stop();
    merge(std::move(temp_writers), components);
}

} // namespace

bool is_multithreaded(const GenomeCallingComponents& components)
{
    return !components.num_threads() || *components.num_threads() > 1;
}

void run_calling(GenomeCallingComponents& components)
{
    if (is_multithreaded(components)) {
        if (DEBUG_MODE) {
            logging::WarningLogger warn_log {};
            warn_log << "Running in parallel mode can make debug log difficult to interpret";
        }
        run_octopus_multi_threaded(components);
    } else {
        run_octopus_single_threaded(components);
    }
}

void destroy(VcfWriter& writer)
{
    VcfWriter tmp {};
    swap(writer, tmp);
}

void log_filtering_info(const GenomeCallingComponents& components)
{
    logging::InfoLogger log {};
    log << "Starting Call Set Refinement (CSR) filtering";
}

std::vector<GenomicRegion> extract_call_regions(VcfReader& vcf)
{
    std::deque<GenomicRegion> regions {};
    auto p = vcf.iterate(VcfReader::UnpackPolicy::sites);
    std::transform(std::move(p.first), std::move(p.second), std::back_inserter(regions),
                   [] (const VcfRecord& record) { return mapped_region(record); });
    return {std::make_move_iterator(std::begin(regions)), std::make_move_iterator(std::end(regions))};
}

std::vector<GenomicRegion> extract_call_regions(boost::filesystem::path vcf_path)
{
    VcfReader tmp {std::move(vcf_path)};
    return extract_call_regions(tmp);
}

template <typename Map>
std::size_t sum_mapped_container_size(const Map& map)
{
    return std::accumulate(std::cbegin(map), std::cend(map), std::size_t {0},
                           [] (auto curr, const auto& p) noexcept { return curr + p.second.size(); });
}

std::vector<GenomicRegion> flatten(const InputRegionMap& regions)
{
    std::vector<GenomicRegion> result {};
    result.reserve(sum_mapped_container_size(regions));
    for (const auto& p : regions) {
        std::copy(std::cbegin(p.second), std::cend(p.second), std::back_inserter(result));
    }
    return result;
}

bool use_unfiltered_call_region_hints_for_filtering(const GenomeCallingComponents& components)
{
    // TODO: no need to do this if reads won't be used, or if likely hints won't help because the calls
    // are very dense.
    return true;
}

HaplotypeLikelihoodModel make_filtering_haplotype_likelihood_model(const GenomeCallingComponents& components)
{
    return components.realignment_haplotype_likelihood_model();
}

Pedigree get_pedigree(const GenomeCallingComponents& components)
{
    auto result = components.pedigree();
    if (!result) {
        result = Pedigree {components.samples().size()};
        for (const auto& sample : components.samples()) {
            result->add_founder({sample});
        }
    }
    return *result;
}

void run_csr(GenomeCallingComponents& components)
{
    if (apply_csr(components)) {
        log_filtering_info(components);
        ProgressMeter progress {components.search_regions()};
        const auto& filter_factory = components.call_filter_factory();
        const auto& filter_read_pipe = components.filter_read_pipe();
        boost::optional<boost::filesystem::path> input_path {};
        if (components.filter_request()) {
            input_path = components.filter_request();
        } else {
            input_path = components.output().path();
        }
        assert(input_path); // cannot be stdout
        BufferedReadPipe::Config buffer_config {components.read_buffer_size()};
        buffer_config.fetch_expansion = 100;
        buffer_config.max_hint_gap = 5'000;
        BufferedReadPipe buffered_rp {filter_read_pipe, buffer_config};
        if (use_unfiltered_call_region_hints_for_filtering(components)) {
            buffered_rp.hint(extract_call_regions(*input_path));
        } else {
            buffered_rp.hint(flatten(components.search_regions()));
        }
        const VcfReader in {std::move(*input_path), components.reference()};
        const auto filter = filter_factory.make(components.reference(), std::move(buffered_rp), in.fetch_header(),
                                                components.ploidies(),
                                                make_filtering_haplotype_likelihood_model(components),
                                                get_pedigree(components),
                                                progress, components.num_threads());
        assert(filter);
        VcfWriter& out {*components.filtered_output()};
        boost::optional<VcfHeader> template_header {};
        namespace fs = boost::filesystem;
        if (components.filter_request() && !components.output().is_open() && components.output().path() && fs::exists(*components.output().path())) {
            template_header = read_vcf_header(*components.output().path());
        }
        filter->filter(in, out, template_header);
        out.close();
    }
}

void log_run_start(const GenomeCallingComponents& components, const UserCommandInfo& info)
{
    static auto debug_log = get_debug_log();
    log_startup_info(components);
    if (debug_log) {
        stream(*debug_log) << "Command line: " << info.command;
        print_input_regions(stream(*debug_log), components.search_regions());
        const auto reads_profile = components.reads_profile();
        if (reads_profile) stream(*debug_log) << "Reads profile:\n" << *reads_profile;
    }
}

class CallingBug : public ProgramError
{
    std::string do_where() const override { return "run_octopus"; }
    std::string do_why() const override
    {
        if (what_) {
            return "Encountered an exception during calling '" + *what_ + "'. This means there is a bug"
                                                                          " and your results are untrustworthy.";
        } else {
            return "Encountered an unknown error during calling. This means there is a bug"
                   " and your results are untrustworthy.";
        }
    }
    
    boost::optional<std::string> what_;
public:
    CallingBug() = default;
    CallingBug(const std::exception& e) : what_ {e.what()} {}
};

void run_variant_calling(GenomeCallingComponents& components, UserCommandInfo info)
{
    static auto debug_log = get_debug_log();
    log_run_start(components, info);
    write_caller_output_header(components, info);
    const auto start = std::chrono::system_clock::now();
    try {
        if (!components.filter_request()) {
            run_calling(components);
        }
    } catch (const Error& e) {
        try {
            if (debug_log) *debug_log << "Encountered an error whilst calling, attempting to cleanup";
            cleanup(components);
        } catch (...) {}
        throw;
    } catch (const std::exception& e) {
        try {
            if (debug_log) *debug_log << "Encountered an error whilst calling, attempting to cleanup";
            cleanup(components);
        } catch (...) {}
        throw CallingBug {e};
    } catch (...) {
        try {
            if (debug_log) *debug_log << "Encountered an error whilst calling, attempting to cleanup";
            cleanup(components);
        } catch (...) {}
        throw CallingBug {};
    }
    components.output().close();
    try {
        run_csr(components);
    } catch (const Error& e) {
        try {
            if (debug_log) *debug_log << "Encountered an error whilst filtering, attempting to cleanup";
            cleanup(components);
        } catch (...) {}
        throw;
    } catch (const std::exception& e) {
        try {
            if (debug_log) *debug_log << "Encountered an error whilst filtering, attempting to cleanup";
            cleanup(components);
        } catch (...) {}
        throw CallingBug {e};
    } catch (...) {
        try {
            if (debug_log) *debug_log << "Encountered an error whilst filtering, attempting to cleanup";
            cleanup(components);
        } catch (...) {}
        throw CallingBug {};
    }
    const auto end = std::chrono::system_clock::now();
    log_finish_info(components, {start, end});
}

bool is_bam_realignment_requested(const GenomeCallingComponents& components)
{
    return static_cast<bool>(components.bamout());
}

bool is_stdout_final_output(const GenomeCallingComponents& components)
{
    return (components.filtered_output() && !components.filtered_output()->path()) || !components.output().path();
}

bool check_bam_realign(const GenomeCallingComponents& components)
{
    logging::WarningLogger warn_log {};
    if (!components.read_manager().all_readers_have_one_sample()) {
        warn_log << "BAM realignment currently only supported for single sample BAMs";
        return false;
    }
    if (is_stdout_final_output(components)) {
        warn_log << "BAM realignment is not supported for stdout calling";
        return false;
    }
    return true;
}

auto get_bam_realignment_vcf(const GenomeCallingComponents& components)
{
    if (components.filtered_output()) {
        return *components.filtered_output()->path();
    } else {
        return *components.output().path();
    }
}

bool is_sam_type(const boost::filesystem::path& path)
{
    const auto type = path.extension().string();
    return type == ".bam" || type == ".sam" || type == ".cram";
}

bool is_called_ploidy_known(const GenomeCallingComponents& components)
{
    const auto contigs = components.contigs();
    return std::all_of(std::cbegin(contigs), std::cend(contigs), [&] (const auto& contig) {
        const auto caller = components.caller_factory().make(components.contigs().front());
        return caller->min_callable_ploidy() == caller->max_callable_ploidy();
    });
}

auto get_max_called_ploidy(VcfReader& vcf, const std::vector<VcfRecord::SampleName>& samples)
{
    unsigned result {0};
    for (auto p = vcf.iterate(); p.first != p.second; ++p.first) {
        for (const auto& sample : samples) {
            result = std::max(p.first->ploidy(sample), result);
        }
    }
    return result;
}

auto get_max_called_ploidy(VcfReader& vcf)
{
    return get_max_called_ploidy(vcf, vcf.fetch_header().samples());
}

auto get_max_called_ploidy(const boost::filesystem::path& output_vcf)
{
    VcfReader vcf {output_vcf};
    return get_max_called_ploidy(vcf);
}

auto get_bam_sampltes(const boost::filesystem::path& bam_path)
{
    io::ReadReader bam {bam_path};
    return bam.extract_samples();
}

auto get_max_called_ploidy(const boost::filesystem::path& output_vcf, const boost::filesystem::path& in_bam)
{
    auto bam_samples = get_bam_sampltes(in_bam);
    std::sort(std::begin(bam_samples), std::end(bam_samples));
    VcfReader vcf {output_vcf};
    auto vcf_samples = vcf.fetch_header().samples();
    std::sort(std::begin(vcf_samples), std::end(vcf_samples));
    std::vector<VcfRecord::SampleName> usable_samples {};
    usable_samples.reserve(std::min(bam_samples.size(), vcf_samples.size()));
    std::set_intersection(std::cbegin(bam_samples), std::cend(bam_samples),
                          std::cbegin(vcf_samples), std::cend(vcf_samples),
                          std::back_inserter(usable_samples));
    return get_max_called_ploidy(vcf, usable_samples);
}

auto get_max_ploidy(const GenomeCallingComponents& components)
{
    if (is_called_ploidy_known(components)) {
        return get_max_ploidy(components.samples(), components.contigs(), components.ploidies());
    } else {
        assert(get_final_output_path(components));
        return get_max_called_ploidy(*get_final_output_path(components));
    }
}

auto get_haplotype_bam_paths(const boost::filesystem::path& prefix, const unsigned max_ploidy)
{
    std::vector<boost::filesystem::path> result {};
    result.reserve(max_ploidy + 1); // + 1 for unassigned reads
    for (unsigned i {1}; i <= max_ploidy + 1; ++i) {
        result.emplace_back(prefix.string() + "_" + std::to_string(i) + ".bam");
    }
    return result;
}

void run_bam_realign(GenomeCallingComponents& components)
{
    if (is_bam_realignment_requested(components)) {
        if (check_bam_realign(components)) {
            components.read_manager().close();
            if (components.read_manager().paths().size() == 1) {
                realign(components.read_manager().paths().front(), get_bam_realignment_vcf(components),
                        *components.bamout(), components.reference(), components.bamout_config());
            } else {
                namespace fs = boost::filesystem;
                const auto bamout_directory = *components.bamout();
                if (fs::exists(bamout_directory)) {
                    if (!fs::is_directory(bamout_directory)) {
                        logging::ErrorLogger error_log {};
                        stream(error_log) << "The given evidence bam directory " << bamout_directory << " is not a directory";
                        return;
                    }
                } else {
                    if (!fs::create_directory(bamout_directory)) {
                        logging::ErrorLogger error_log {};
                        stream(error_log) << "Failed to create temporary directory " << bamout_directory << " - check permissions";
                        return;
                    }
                }
                for (const auto& bamin_path : components.read_manager().paths()) {
                    auto bamout_path = bamout_directory;
                    bamout_path /= bamin_path.filename();
                    if (bamin_path != bamout_path) {
                        realign(bamin_path, get_bam_realignment_vcf(components), bamout_path, components.reference(), components.bamout_config());
                    } else {
                        logging::WarningLogger warn_log {};
                        stream(warn_log) << "Cannot make evidence bam " << bamout_path << " as it is an input bam";
                    }
                }
            }
        }
    }
}

void run_data_profiler(GenomeCallingComponents& components)
{
    const auto data_profile_csv_path = components.data_profile();
    if (data_profile_csv_path) {
        logging::InfoLogger info_log {};
        VcfWriter& final_output {get_final_output(components)};
        const auto final_output_path = final_output.path();
        if (final_output_path) {
            info_log << "Starting indel profiler";
            final_output.close();
            if (is_indexable(*final_output_path)) index_vcf(*final_output_path);
            auto config = components.profiler_config();
            const auto profile = profile_indels(components.read_pipe(), *final_output_path, components.reference(), components.search_regions(), std::move(config));
            std::ofstream profile_file {data_profile_csv_path->string()};
            profile_file << profile;
            stream(info_log) << "Indel profile written to " << *data_profile_csv_path;
        } else {
            info_log << "Did not run indel profiler as calls not written to file";
        }
    }
}

void run_post_calling_requests(GenomeCallingComponents& components)
{
    run_data_profiler(components);
    run_bam_realign(components);
}

void run_octopus(GenomeCallingComponents& components, UserCommandInfo info)
{
    run_variant_calling(components, std::move(info));
    run_post_calling_requests(components);
    cleanup(components);
}

} // namespace octopus
