// Copyright (c) 2016 Daniel Cooke
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

#include "config/common.hpp"
#include "basics/genomic_region.hpp"
#include "concepts/mappable.hpp"
#include "containers/mappable_flat_multi_set.hpp"
#include "containers/mappable_map.hpp"
#include "io/reference/reference_genome.hpp"
#include "io/read/read_manager.hpp"
#include "readpipe/read_pipe_fwd.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/read_stats.hpp"
#include "config/octopus_vcf.hpp"
#include "core/callers/caller_factory.hpp"
#include "core/callers/caller.hpp"
#include "utils/maths.hpp"
#include "logging/progress_meter.hpp"
#include "logging/logging.hpp"
#include "core/callers/utils/vcf_header_factory.hpp"
#include "io/variant/vcf.hpp"
#include "utils/timing.hpp"

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

using CallTypeSet = std::set<std::type_index>;

VcfHeader make_vcf_header(const std::vector<SampleName>& samples,
                          const std::vector<GenomicRegion::ContigName>& contigs,
                          const ReferenceGenome& reference,
                          const CallTypeSet& call_types,
                          const std::string& command)
{
    auto builder = vcf::make_header_template().set_samples(samples);
    for (const auto& contig : contigs) {
        builder.add_contig(contig, {{"length", std::to_string(reference.contig_size(contig))}});
    }
    builder.add_basic_field("reference", reference.name());
    builder.add_structured_field("octopus", {{"command", '"' + command + '"'}});
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
                          const std::string& command)
{
    return make_vcf_header(samples, std::vector<GenomicRegion::ContigName> {contig},
                           reference, call_types, command);
}

bool has_reads(const GenomicRegion& region, ContigCallingComponents& components)
{
    return components.read_manager.get().has_reads(components.samples.get(), region);
}

auto get_call_types(const GenomeCallingComponents& components, const std::vector<ContigName>& contigs)
{
    CallTypeSet result {};
    
    for (const auto& contig : components.contigs()) {
        const auto tmp_caller = components.caller_factory().make(contig);
        auto caller_call_types = tmp_caller->call_types();
        result.insert(std::begin(caller_call_types), std::end(caller_call_types));
    }
    
    return result;
}

void write_caller_output_header(GenomeCallingComponents& components,
                                const std::string& command)
{
    const auto call_types = get_call_types(components, components.contigs());
    if (components.sites_only()) {
        components.output() << make_vcf_header({}, components.contigs(), components.reference(),
                                               call_types, command);
    } else {
        components.output() << make_vcf_header(components.samples(), components.contigs(),
                                               components.reference(), call_types, command);
    }
}

std::string get_caller_name(const GenomeCallingComponents& components)
{
    assert(!components.contigs().empty());
    const auto& test_contig = components.contigs().front();
    const auto test_caller = components.caller_factory().make(test_contig);
    return test_caller->name(); // There can only be one caller type per run
}

void log_startup_info(const GenomeCallingComponents& components)
{
    logging::InfoLogger log {};
    stream(log) << "Invoked calling model: " << get_caller_name(components);
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
    auto sl = stream(log);
    sl << "Writing calls to ";
    const auto output_path = components.output().path();
    if (output_path) {
        sl << *output_path;
    } else {
        sl << "stdout";
    }
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
                            const GenomicRegion& remaining_call_region,
                            boost::optional<GenomicRegion::Size> min_size = boost::none)
{
    if (is_empty(remaining_call_region)) {
        return remaining_call_region;
    }
    const auto max_window = find_max_window(components, remaining_call_region);
    if (ends_before(remaining_call_region, max_window)) {
        return remaining_call_region;
    }
    if (min_size && size(max_window) < *min_size) {
        if (size(remaining_call_region) < *min_size) {
            return remaining_call_region;
        }
        return expand_rhs(head_region(max_window), *min_size);
    }
    return max_window;
}

auto propose_call_subregion(const ContigCallingComponents& components,
                            const GenomicRegion& current_subregion,
                            const GenomicRegion& input_region,
                            boost::optional<GenomicRegion::Size> min_size = boost::none)
{
    assert(contains(input_region, current_subregion));
    return propose_call_subregion(components, right_overhang_region(input_region, current_subregion), min_size);
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
    for (const auto& contig : components.contigs()) {
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
    boost::filesystem::path file_name {contig + "_" + begin + "-" + end + "_temp.bcf"};
    path /= file_name;
    const auto call_types = get_call_types(components, {region.contig_name()});
    auto header = make_vcf_header(components.samples(), contig, components.reference(), call_types,
                                  "octopus-internal");
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
    result.reserve(components.contigs().size());
    for (const auto& contig : components.contigs()) {
        result.emplace(contig, create_unique_temp_output_file(contig, components));
    }
    return result;
}

struct Task : public Mappable<Task>
{
    enum class ExecutionPolicy { seq, par, par_vec }; // To match Parallelism TS
    
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

TaskQueue divide_work_into_tasks(const ContigCallingComponents& components,
                                 const Task::ExecutionPolicy policy)
{
    TaskQueue result {};
    
    if (components.regions.empty()) return result;
    
    static constexpr GenomicRegion::Size minTaskSize {1000};
    
    for (const auto& region : components.regions) {
        auto subregion = propose_call_subregion(components, region, minTaskSize);
        do {
            result.emplace(subregion, policy);
            subregion = propose_call_subregion(components, subregion, region, minTaskSize);
        } while (!is_after(subregion, region));
    }
    
    return result;
}

Task::ExecutionPolicy make_execution_policy(const GenomeCallingComponents& components)
{
    if (components.num_threads()) {
        return Task::ExecutionPolicy::seq;
    }
    return Task::ExecutionPolicy::par;
}

struct TaskMakerSyncPacket
{
    std::condition_variable cv;
    std::mutex mutex;
    std::atomic_bool done;
};

auto make_contig_components(const ContigName& contig, GenomeCallingComponents& components,
                            const unsigned num_threads)
{
    ContigCallingComponents result {contig, components};
    result.read_buffer_size /= num_threads;
    return result;
}

void make_remaining_tasks(TaskMap& tasks, std::vector<ContigName> contigs, GenomeCallingComponents& components,
                          const unsigned num_threads, Task::ExecutionPolicy policy, TaskMakerSyncPacket& sync)
{
    assert(!contigs.empty());
    std::unique_lock<std::mutex> lk {sync.mutex, std::defer_lock};
    
    std::for_each(std::make_move_iterator(std::begin(contigs)),
                  std::make_move_iterator(std::prev(std::end(contigs))),
                  [&] (auto&& contig) {
                      auto contig_components = make_contig_components(contig, components, num_threads);
                      auto contig_tasks = divide_work_into_tasks(contig_components, policy);
                      
                      lk.lock();
                      tasks.emplace(std::move(contig), std::move(contig_tasks));
                      lk.unlock();
                      sync.cv.notify_all();
                  });
    auto contig_components = make_contig_components(contigs.back(), components, num_threads);
    auto contig_tasks = divide_work_into_tasks(contig_components, policy);
    
    lk.lock();
    tasks.emplace(std::move(contigs.back()), std::move(contig_tasks));
    sync.done = true;
    lk.unlock();
    sync.cv.notify_all();
}

std::thread make_tasks(TaskMap& tasks, GenomeCallingComponents& components, const unsigned num_threads,
                       TaskMakerSyncPacket& sync)
{
    auto contigs = components.contigs();
    if (contigs.empty()) {
        sync.done = true;
        return std::thread {};
    }
    const auto policy = make_execution_policy(components);
    auto contig_components = make_contig_components(contigs.front(), components, num_threads);
    tasks.emplace(std::move(contigs.front()), divide_work_into_tasks(contig_components, policy));
    if (contigs.size() == 1) {
        sync.done = true;
        return std::thread {};
    }
    contigs.erase(std::cbegin(contigs));
    return std::thread {
        make_remaining_tasks, std::ref(tasks), std::move(contigs), std::ref(components),
        num_threads, policy, std::ref(sync)
    };
}

void log_num_cores(const unsigned num_cores)
{
    auto debug_log = logging::get_debug_log();
    if (debug_log) stream(*debug_log) << "Detected " << num_cores << " system cores";
}

void warn_undetected_cores()
{
    logging::WarningLogger log {};
    log << "Unable to detect the number of system cores,"
    " it may be better to run with a user number if the number of cores is known";
}

unsigned calculate_num_task_threads(const GenomeCallingComponents& components)
{
    if (components.num_threads()) {
        return *components.num_threads();
    }
    // TODO: come up with a better calculation
    const auto num_cores = std::thread::hardware_concurrency();
    if (num_cores > 0) {
        log_num_cores(num_cores);
        return num_cores;
    } else {
        warn_undetected_cores();
        return std::min(components.read_manager().num_files(), 8u);
    }
}

Task pop(TaskMap& tasks)
{
    assert(!tasks.empty());
    const auto it = std::begin(tasks);
    assert(!it->second.empty());
    const auto result = std::move(it->second.front());
    it->second.pop();
    if (it->second.empty()) {
        tasks.erase(it);
    }
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

struct SyncPacket
{
    std::condition_variable cv;
    std::mutex mutex;
    std::atomic_uint count_finished;
};

template<typename R>
bool is_ready(const std::future<R>& f)
{
    return f.valid() && f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
}

auto run(Task task, ContigCallingComponents components, SyncPacket& sync)
{
    static auto debug_log = get_debug_log();
    if (debug_log) stream(*debug_log) << "Spawning task " << task;
    
    return std::async(std::launch::async,
                      [task = std::move(task), components = std::move(components), &sync]
                      () -> CompletedTask {
                          CompletedTask result {task};
                          result.runtime.start = std::chrono::system_clock::now();
                          result.calls = components.caller->call(task.region, components.progress_meter);
                          result.runtime.end = std::chrono::system_clock::now();
                          
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
                        [&rhs_begin] (const auto& call) {
                            return mapped_end(call) > rhs_begin;
                        });
}

auto find_last_rhs_connecting(const GenomicRegion& lhs_region, const std::deque<VcfRecord>& rhs_calls)
{
    const auto lhs_end = mapped_end(lhs_region);
    return std::find_if_not(std::cbegin(rhs_calls), std::cend(rhs_calls),
                            [&lhs_end] (const auto& call) {
                                return mapped_begin(call) < lhs_end;
                            });
}

void resolve_connecting_calls(CompletedTask& lhs, CompletedTask& rhs,
                              const ContigCallingComponentFactory& calling_components)
{
    static auto debug_log = get_debug_log();
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::make_move_iterator;
    if (lhs.calls.empty() || rhs.calls.empty()) return;
    const auto first_lhs_connecting = find_first_lhs_connecting(lhs.calls, encompassing_region(rhs.calls));
    const auto last_rhs_connecting  = find_last_rhs_connecting(encompassing_region(lhs.calls), rhs.calls);
    if (first_lhs_connecting == cend(lhs.calls) && last_rhs_connecting == cbegin(rhs.calls)) {
        return;
    }
    if (debug_log) {
        stream(*debug_log) << "Resolving connecting calls between tasks " << lhs << " & " << rhs;
    }
    std::deque<VcfRecord> merged_calls {};
    std::set_union(first_lhs_connecting, cend(lhs.calls),
                   cbegin(rhs.calls), last_rhs_connecting,
                   std::back_inserter(merged_calls));
    lhs.calls.erase(first_lhs_connecting, cend(lhs.calls));
    rhs.calls.erase(cbegin(rhs.calls), last_rhs_connecting);
    
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
                if (!contains(unresolved_region, encompassing_region(resolved_calls))) {
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
    if (adjacent_tasks.size() < 2) return;
    assert(std::is_sorted(std::cbegin(adjacent_tasks), std::cend(adjacent_tasks)));
    auto lhs = std::begin(adjacent_tasks);
    std::for_each(std::next(lhs), std::end(adjacent_tasks),
                  [&] (auto& rhs) { resolve_connecting_calls(*lhs++, rhs, calling_components); });
}

void write(std::deque<CompletedTask>&& tasks, VcfWriter& temp_vcf)
{
    static auto debug_log = get_debug_log();
    for (auto&& task : tasks) {
        if (debug_log) {
            stream(*debug_log) << "Writing completed task " << task
                               << " that finished in " << duration(task);
        }
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
        
        // Keep the last task buffered to enable connection resolution when the next task finishes
        
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

void extract_remaining_future_tasks(FutureCompletedTasks& futures, std::deque<CompletedTask>& result)
{
    const auto it = std::remove_if(std::begin(futures), std::end(futures),
                                   [] (const auto& f) { return !f.valid(); });
    std::transform(std::begin(futures), it, std::back_inserter(result),
                   [] (auto& fut) { return fut.get(); });
    futures.clear();
    futures.shrink_to_fit();
}

void extract_buffered_tasks(CompletedTaskMap& buffered_tasks, std::deque<CompletedTask>& result)
{
    for (auto& p : buffered_tasks) {
        std::transform(std::make_move_iterator(std::begin(p.second)),
                       std::make_move_iterator(std::end(p.second)),
                       std::back_inserter(result),
                       [] (auto&& p) { return std::move(p.second); });
        p.second.clear();
    }
    buffered_tasks.clear();
}

void sort_ignoring_contig_name(std::deque<CompletedTask>& tasks)
{
    std::sort(std::begin(tasks), std::end(tasks),
              [] (const auto& lhs, const auto& rhs) {
                  return contig_region(lhs) < contig_region(rhs);
              });
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

void write_remaining_tasks(FutureCompletedTasks& futures, CompletedTaskMap& buffered_tasks,
                           TempVcfWriterMap& temp_vcfs,
                           const ContigCallingComponentFactoryMap& calling_components)
{
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
    return writers_to_readers(extract_writers(std::move(vcfs)));
}

void merge(TempVcfWriterMap&& temp_vcf_writers, GenomeCallingComponents& components)
{
    auto temp_readers = extract_as_readers(std::move(temp_vcf_writers));
    merge(temp_readers, components.output(), components.contigs());
}

bool done(const TaskMap& pending_tasks, TaskMakerSyncPacket& sync) noexcept
{
    std::lock_guard<std::mutex> lk {sync.mutex};
    return pending_tasks.empty() && sync.done;
}

void run_octopus_multi_threaded(GenomeCallingComponents& components)
{
    using namespace std::chrono_literals;
    static auto debug_log = get_debug_log();
    
    auto temp_vcfs = make_temp_vcf_writers(components);
    const auto num_task_threads = calculate_num_task_threads(components);
    
    TaskMap pending_tasks {components.contigs()};
    TaskMakerSyncPacket task_maker_sync {};
    std::unique_lock<std::mutex> pending_task_lock {task_maker_sync.mutex, std::defer_lock};
    auto maker_thread = make_tasks(pending_tasks, components, num_task_threads, task_maker_sync);
    if (maker_thread.joinable()) maker_thread.detach();
    
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
    
    SyncPacket task_sync {};
    const auto calling_components = make_contig_calling_component_factory_map(components);
    unsigned num_idle_futures {0};
    
    components.progress_meter().start();
    
    while (!done(pending_tasks, task_maker_sync)) {
        pending_task_lock.lock();
        if (!task_maker_sync.done && pending_tasks.empty()) {
            if (num_idle_futures < futures.size()) {
                // If there are running futures then it's good periodically check to see if
                // any have finished and process them while we wait for the task maker.
                while (pending_tasks.empty() && task_sync.count_finished == 0) {
                    auto now = std::chrono::system_clock::now();
                    task_maker_sync.cv.wait_until(pending_task_lock, now + 5s, [&] () { return !pending_tasks.empty(); });
                }
            } else {
                task_maker_sync.cv.wait(pending_task_lock, [&] () { return !pending_tasks.empty(); });
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
                                temp_vcfs.at(contig), calling_components.at(contig));
                --task_sync.count_finished;
            }
            if (!future.valid()) {
                pending_task_lock.lock();
                if (!pending_tasks.empty()) {
                    auto task = pop(pending_tasks);
                    pending_task_lock.unlock();
                    future = run(task, calling_components.at(contig_name(task))(), task_sync);
                    running_tasks.at(contig_name(task)).push(std::move(task));
                } else {
                    pending_task_lock.unlock();
                    ++num_idle_futures;
                }
            }
        }
        // If there are idle futures then we must have run out of tasks, so we should
        // wait on the task maker condition variable instead.
        if (num_idle_futures == 0 && task_sync.count_finished == 0) {
            std::unique_lock<std::mutex> lk {task_sync.mutex};
            task_sync.cv.wait(lk, [&] () { return task_sync.count_finished > 0; });
        } else {
            if (debug_log) stream(*debug_log) << "There are " << num_idle_futures << " idle futures";
        }
    }
    running_tasks.clear();
    holdbacks.clear(); // holdbacks are just references to buffered tasks
    write_remaining_tasks(futures, buffered_tasks, temp_vcfs, calling_components);
    components.progress_meter().stop();
    merge(std::move(temp_vcfs), components);
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
    components.output().close();
}

auto make_filter_read_pipe(const GenomeCallingComponents& components)
{
    using std::make_unique;
    using namespace readpipe;
    
    ReadTransformer transformer {};
    transformer.register_transform(MaskSoftClipped {});
    
    using ReadFilterer = ReadPipe::ReadFilterer;
    ReadFilterer filterer {};
    filterer.add(make_unique<HasValidBaseQualities>());
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

auto add_identifier(const boost::filesystem::path& base, const std::string& identifier)
{
    const auto old_stem  = base.stem();
    const auto extension = base.extension();
    
    boost::filesystem::path new_stem;
    
    if (extension.string() == ".gz") {
        new_stem = old_stem.stem().string() + "." + identifier
        + old_stem.extension().string() + extension.string();
    } else {
        new_stem = old_stem.string() + "." + identifier + extension.string();
    }
    
    return base.parent_path() / new_stem;
}

void run_filtering(const GenomeCallingComponents& components)
{
    auto filter = components.call_filter();
    if (filter) {
        // TODO
    }
}

auto get_legacy_path(const boost::filesystem::path& native)
{
    return add_identifier(native, "legacy");
}

void run_reporting(const GenomeCallingComponents& components)
{
    if (components.legacy()) {
        auto output_path = components.output().path();
        if (output_path) {
            const VcfReader native {*output_path};
            VcfWriter legacy {get_legacy_path(native.path())};
            convert_to_legacy(native, legacy);
        }
    }
}

void log_run_start(const GenomeCallingComponents& components)
{
    static auto debug_log = get_debug_log();
    log_startup_info(components);
    if (debug_log) print_input_regions(stream(*debug_log), components.search_regions());
}

void run_octopus(GenomeCallingComponents& components, std::string command)
{
    static auto debug_log = get_debug_log();
    logging::InfoLogger info_log {};
    using utils::TimeInterval;
    
    const auto start = std::chrono::system_clock::now();
    try {
        log_run_start(components);
        write_caller_output_header(components, command);
        run_calling(components);
        run_filtering(components);
        run_reporting(components);
    } catch (const std::exception& e) {
        try {
            if (debug_log) *debug_log << "Encountered an error, attempting to cleanup";
            cleanup(components);
        } catch (...) {}
        throw;
    }
    const auto end = std::chrono::system_clock::now();
    const auto search_size = sum_region_sizes(components.search_regions());
    stream(info_log) << "Finished calling "
                     << utils::format_with_commas(search_size) << "bp, total runtime "
                     << TimeInterval {start, end};
    cleanup(components);
}
} // namespace octopus


