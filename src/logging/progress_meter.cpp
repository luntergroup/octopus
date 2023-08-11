// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "progress_meter.hpp"

#include <sstream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <utility>
#include <initializer_list>
#include <cassert>

#include "utils/mappable_algorithms.hpp"
#include "containers/mappable_map.hpp"
#include "utils/timing.hpp"
#include "utils/string_utils.hpp"
#include "utils/maths.hpp"

namespace octopus {

using utils::TimeInterval;

namespace {

template <typename T>
unsigned num_digits(const T x)
{
    return static_cast<unsigned>(std::to_string(x).size());
}

template <typename T>
unsigned max_str_length(const T& p)
{
    auto result = static_cast<unsigned>(p.first.size());
    if (p.second.empty()) return result;
    return result + num_digits(p.second.rightmost().end());
}

unsigned max_position_str_length(const InputRegionMap& input_regions)
{
    assert(!input_regions.empty());
    unsigned result {0};
    for (const auto& p : input_regions) {
        const auto curr = max_str_length(p);
        if (result < curr) result = curr;
    }
    return result;
}

unsigned max_position_str_length(const GenomicRegion& region)
{
    return static_cast<unsigned>(region.contig_name().size()) + num_digits(region.end());
}

template <typename T>
unsigned calculate_position_tab_length(const T& region)
{
    return std::max(18u, max_position_str_length(region));
}

} // namespace

ProgressMeter::ProgressMeter(InputRegionMap regions)
: target_regions_ {std::move(regions)}
, completed_regions_ {}
, num_bp_to_search_ {}
, num_bp_completed_ {0}
, percent_until_tick_ {max_tick_size_}
, percent_at_last_tick_ {0}
, start_ {std::chrono::system_clock::now()}
, last_tick_ {start_}
, tick_durations_ {}
, done_ {false}
, position_tab_length_ {}
, block_compute_times_ {}
, log_ {}
, buffer_size_ {}
, buffer_ {}
{
    for (auto& p : target_regions_) {
        auto covered_regions = extract_covered_regions(p.second);
        p.second.clear();
        p.second.insert(std::make_move_iterator(std::begin(covered_regions)),
                        std::make_move_iterator(std::end(covered_regions)));
        p.second.shrink_to_fit();
    }
    num_bp_to_search_ = sum_region_sizes(target_regions_);
    if (!target_regions_.empty()) {
        position_tab_length_ = calculate_position_tab_length(target_regions_);
    }
}

ProgressMeter::ProgressMeter(GenomicRegion region)
: ProgressMeter {InputRegionMap {std::make_pair(region.contig_name(), InputRegionMap::mapped_type {region})}}
{}

ProgressMeter::ProgressMeter(ProgressMeter&& other)
{
    std::lock_guard<std::mutex> lock {other.mutex_};
    target_regions_       = std::move(other.target_regions_);
    completed_regions_    = std::move(other.completed_regions_);
    num_bp_to_search_     = std::move(other.num_bp_to_search_);
    num_bp_completed_     = std::move(other.num_bp_completed_);
    max_tick_size_        = std::move(other.max_tick_size_);
    curr_tick_size_       = std::move(other.curr_tick_size_);
    percent_until_tick_   = std::move(other.percent_until_tick_);
    percent_at_last_tick_ = std::move(other.percent_at_last_tick_);
    start_                = std::move(other.start_);
    last_tick_            = std::move(other.last_tick_);
    tick_durations_       = std::move(other.tick_durations_);
    done_                 = std::move(other.done_);
    position_tab_length_  = std::move(other.position_tab_length_);
    block_compute_times_  = std::move(other.block_compute_times_);
    log_                  = std::move(other.log_);
    buffer_size_          = std::move(other.buffer_size_);
    buffer_               = std::move(other.buffer_);
}

ProgressMeter& ProgressMeter::operator=(ProgressMeter&& other)
{
    if (this != &other) {
        std::unique_lock<std::mutex> lock_lhs {mutex_, std::defer_lock}, lock_rhs {other.mutex_, std::defer_lock};
        std::lock(lock_lhs, lock_rhs);
        target_regions_       = std::move(other.target_regions_);
        completed_regions_    = std::move(other.completed_regions_);
        num_bp_to_search_     = std::move(other.num_bp_to_search_);
        num_bp_completed_     = std::move(other.num_bp_completed_);
        max_tick_size_        = std::move(other.max_tick_size_);
        curr_tick_size_       = std::move(other.curr_tick_size_);
        percent_until_tick_   = std::move(other.percent_until_tick_);
        percent_at_last_tick_ = std::move(other.percent_at_last_tick_);
        start_                = std::move(other.start_);
        last_tick_            = std::move(other.last_tick_);
        tick_durations_       = std::move(other.tick_durations_);
        done_                 = std::move(other.done_);
        position_tab_length_  = std::move(other.position_tab_length_);
        block_compute_times_  = std::move(other.block_compute_times_);
        log_                  = std::move(other.log_);
        buffer_size_          = std::move(other.buffer_size_);
        buffer_               = std::move(other.buffer_);
    }
    return *this;
}

namespace {

double percent_completed(const std::size_t num_bp_completed, const std::size_t num_bp_to_search)
{
    return 100 * static_cast<double>(num_bp_completed) / num_bp_to_search;
}

auto percent_completed_str(const std::size_t num_bp_completed, const std::size_t num_bp_to_search)
{
    return utils::to_string(percent_completed(num_bp_completed, num_bp_to_search), 1) + '%';
}

std::string to_string(const TimeInterval& duration)
{
    std::ostringstream ss{};
    ss << duration;
    return ss.str();
}

template <typename ForwardIt>
auto mean_duration(ForwardIt first, ForwardIt last)
{
    return maths::mean(first, last, [](const auto& d) { return d.count(); });
}

template <typename Container>
auto mean_duration(const Container& durations)
{
    return mean_duration(std::cbegin(durations), std::cend(durations));
}

template <typename ForwardIt>
auto stdev_duration(ForwardIt first, ForwardIt last)
{
    return maths::stdev(first, last, [](const auto& d) { return d.count(); });
}

template <typename Container>
auto stdev_duration(const Container& durations)
{
    return stdev_duration(std::cbegin(durations), std::cend(durations));
}

template <typename Container>
bool all_equal(const Container& c)
{
    return std::adjacent_find(std::cbegin(c), std::cend(c), std::not_equal_to<void> {}) == std::cend(c);
}

template <typename Container>
void remove_outliers(Container& durations)
{
    if (durations.size() < 2 || all_equal(durations)) return;
    auto min_duration_itr = std::min_element(std::begin(durations), std::end(durations));
    if (min_duration_itr == std::begin(durations)) {
        min_duration_itr = std::remove(min_duration_itr, std::end(durations), *min_duration_itr);
    } else {
        min_duration_itr = std::end(durations);
    }
    const auto mean = mean_duration(std::begin(durations), min_duration_itr);
    const auto stdev = stdev_duration(std::begin(durations), min_duration_itr);
    const auto min = std::max(0.0, mean - (2 * stdev));
    const auto max = mean + (2 * stdev);
    min_duration_itr = std::remove_if(min_duration_itr, std::end(durations),
                                      [=](const auto& duration) {
                                          return duration.count() < min || duration.count() > max;
                                      });
    durations.erase(min_duration_itr, std::end(durations));
}

template <typename Container>
auto estimate_ttc(const std::chrono::system_clock::time_point now,
                  const Container& durations,
                  const std::size_t num_remaining_blocks)
{
    const auto mean_block_duration = mean_duration(durations);
    const auto estimated_remaining_ticks = static_cast<std::size_t>(num_remaining_blocks * mean_block_duration);
    const std::chrono::milliseconds estimated_remaining_duration{estimated_remaining_ticks};
    return TimeInterval {now, now + estimated_remaining_duration};
}

} // namespace

ProgressMeter::~ProgressMeter()
{
    if (!done_ && !target_regions_.empty() && num_bp_completed_ > 0) {
        const TimeInterval duration {start_, std::chrono::system_clock::now()};
        const auto time_taken = to_string(duration);
        stream(log_) << std::string(position_tab_length_ - 4, ' ') << "-"
                     << completed_pad("100%", position_tab_length_ - 3) << "100%"
                     << time_taken_pad(time_taken) << time_taken
                     << ttc_pad("-") << "-";
    }
}

void ProgressMeter::set_max_tick_size(double percent)
{
    block_compute_times_.clear(); // TODO: can we use old block times to estimate new block times?
    max_tick_size_ = percent;
    percent_until_tick_  = std::min(percent, percent_until_tick_);
    curr_tick_size_ = std::min(max_tick_size_, curr_tick_size_);
}

void ProgressMeter::set_buffer_size(const std::size_t n)
{
    buffer_size_ = n;
}

void ProgressMeter::start()
{
    if (!target_regions_.empty()) {
        completed_regions_.reserve(target_regions_.size());
        write_header();
    }
    start_ = std::chrono::system_clock::now();
    last_tick_ = start_;
}

void ProgressMeter::resume()
{
    // TODO
}

void ProgressMeter::pause()
{
    // TODO
}

void ProgressMeter::stop()
{
    spill_buffer();
    if (!done_ && !target_regions_.empty()) {
        const TimeInterval duration {start_, std::chrono::system_clock::now()};
        const auto time_taken = to_string(duration);
        stream(log_) << std::string(position_tab_length_ - 4, ' ') << "-"
                     << completed_pad("100%", position_tab_length_ - 3) << "100%"
                     << time_taken_pad(time_taken) << time_taken
                     << ttc_pad("-") << "-";
    }
    done_ = true;
}

void ProgressMeter::reset()
{
    if (!done_) stop();
    completed_regions_.clear();
    num_bp_to_search_ = sum_region_sizes(target_regions_);
    num_bp_completed_ = 0;
    curr_tick_size_ = max_tick_size_;
    percent_until_tick_ = max_tick_size_;
    percent_at_last_tick_ = 0;
    start_ = std::chrono::system_clock::now();
    last_tick_ = start_;
    tick_durations_.clear();
    tick_durations_.shrink_to_fit();
    done_ = false;
    block_compute_times_.clear();
    block_compute_times_.shrink_to_fit();
    buffer_.clear();
    buffer_.shrink_to_fit();
}

void ProgressMeter::log_completed(const GenomicRegion& region)
{
    std::lock_guard<std::mutex> lock {mutex_};
    if (buffer_size_ > 0) {
        buffer(region);
    } else {
        do_log_completed(region);
    }
}

void ProgressMeter::log_completed(const GenomicRegion::ContigName& contig)
{
    if (target_regions_.count(contig) == 1) {
        const auto& contig_regions = target_regions_.at(contig);
        if (!contig_regions.empty()) {
            const auto contig_region = encompassing_region(contig_regions.front(), contig_regions.back());
            log_completed(contig_region);
        }
    }
}

// private methods

void ProgressMeter::buffer(const GenomicRegion& region)
{
    if (!buffer_.empty() && (are_adjacent(buffer_.back(), region) || overlaps(buffer_.back(), region))) {
        buffer_.back() = encompassing_region(buffer_.back(), region);
    } else if (buffer_.size() == buffer_size_) {
        spill_buffer();
        buffer_.assign({region});
    } else {
        buffer_.push_back(region);
    }
}

void ProgressMeter::spill_buffer()
{
    for (const auto& region : buffer_) {
        do_log_completed(region);
    }
}

void ProgressMeter::do_log_completed(const GenomicRegion& region)
{
    const auto new_bp_processed = merge(region);
    const auto new_percent_done = percent_completed(new_bp_processed, num_bp_to_search_);
    num_bp_completed_ += new_bp_processed;
    percent_until_tick_ -= new_percent_done;
    if (percent_until_tick_ <= 0) output_log(region);
}

ProgressMeter::RegionSizeType ProgressMeter::merge(const GenomicRegion& region)
{
    RegionSizeType result {0};
    const auto target_itr = target_regions_.find(region.contig_name());
    if (target_itr != std::cend(target_regions_)) {
        const auto overlapped_targets = overlap_range(target_itr->second, region);
        for (const auto& target : overlapped_targets) {
            const auto completed_target = overlapped_region(target, region)->contig_region();
            if (completed_regions_.count(region.contig_name()) == 0) {
                completed_regions_.emplace(std::piecewise_construct,
                                           std::forward_as_tuple(region.contig_name()),
                                           std::forward_as_tuple(std::initializer_list<ContigRegion>({completed_target})));
                result += size(completed_target);
            } else {
                auto& completed_regions = completed_regions_.at(region.contig_name());
                if (completed_regions.has_overlapped(completed_target)) {
                    auto overlapped = completed_regions.overlap_range(completed_target);
                    assert(size(overlapped) >= 1);
                    if (size(overlapped) == 1) {
                        const auto& interactor = overlapped.front();
                        if (!contains(interactor, completed_target)) {
                            auto new_region = completed_target;
                            if (begins_before(new_region, interactor)) {
                                result += left_overhang_size(new_region, interactor);
                            } else if (begins_before(interactor, new_region)) {
                                new_region = encompassing_region(overlapped.front(), new_region);
                            }
                            if (ends_before(interactor, new_region)) {
                                result += right_overhang_size(new_region, interactor);
                            } else if (ends_before(new_region, interactor)) {
                                new_region = encompassing_region(interactor, new_region);
                            }
                            completed_regions.erase(interactor);
                            completed_regions.insert(std::move(new_region));
                        }
                    } else {
                        auto new_region = completed_target;
                        result += size(new_region);
                        result -= overlap_size(overlapped.front(), completed_target);
                        if (begins_before(overlapped.front(), new_region)) {
                            new_region = encompassing_region(overlapped.front(), new_region);
                        }
                        overlapped.advance_begin(1);
                        result -= overlap_size(overlapped.back(), new_region);
                        if (ends_before(new_region, overlapped.back())) {
                            new_region = encompassing_region(new_region, overlapped.back());
                        }
                        overlapped.advance_end(-1);
                        for (const auto& r : overlapped) result -= size(r);
                        completed_regions.erase_overlapped(completed_target);
                        completed_regions.insert(std::move(new_region));
                    }
                } else {
                    completed_regions.insert(completed_target);
                    result += size(completed_target);
                }
            }
        }
    }
    return result;
}

void ProgressMeter::write_header()
{
    const std::string pos_tab_bar(position_tab_length_, '-');
    assert(position_tab_length_ >= 8);
    const auto num_position_tab_whitespaces = position_tab_length_ - 8;
    const std::string pos_tab_lhs_pad((position_tab_length_ - 8) / 2, ' ');
    auto pos_tab_rhs_pad = pos_tab_lhs_pad;
    if (num_position_tab_whitespaces % 2) {
        pos_tab_rhs_pad += ' ';
    }
    stream(log_) << pos_tab_bar << "------------------------------------------------------";
    stream(log_) << pos_tab_lhs_pad << "current " << pos_tab_rhs_pad
    << "|                   |     time      |     estimated   ";
    stream(log_) << pos_tab_lhs_pad << "position" << pos_tab_rhs_pad
    << "|     completed     |     taken     |     ttc         ";
    stream(log_) << pos_tab_bar << "------------------------------------------------------";
}

void ProgressMeter::output_log(const GenomicRegion& region)
{
    const auto percent_done = percent_completed(num_bp_completed_, num_bp_to_search_);
    const auto now = std::chrono::system_clock::now();
    const TimeInterval duration {start_, now};
    const auto time_taken = to_string(duration);
    if (percent_done >= 100) return;
    const auto percent_since_last_tick = percent_done - percent_at_last_tick_;
    const auto num_blocks_completed = static_cast<std::size_t>(std::floor(percent_since_last_tick / min_tick_size_));
    std::string ttc {"-"};
    if (num_blocks_completed > 0) {
        const auto tick_duration = std::chrono::duration_cast<DurationUnits>(now - last_tick_);
        const DurationUnits duration_per_block {tick_duration.count() / num_blocks_completed};
        std::fill_n(std::back_inserter(block_compute_times_), num_blocks_completed, duration_per_block);
        const auto num_remaining_blocks = static_cast<std::size_t>((100.0 - percent_done) / min_tick_size_);
        remove_outliers(block_compute_times_);
        ttc = to_string(estimate_ttc(now, block_compute_times_, num_remaining_blocks));
        assert(!ttc.empty());
        if (ttc.front() == '0') {
            ttc = "-";
        }
    }
    const auto percent_completed = percent_completed_str(num_bp_completed_, num_bp_to_search_);
    const auto position_tick = region.contig_name() + ":" + std::to_string(region.end());
    const auto position_str = position_pad(region) + position_tick;
    stream(log_) << position_str
                 << completed_pad(percent_completed, position_str.size()) << percent_completed
                 << time_taken_pad(time_taken) << time_taken
                 << ttc_pad(ttc) << ttc;
    tick_durations_.emplace_back(now - last_tick_);
    last_tick_            = now;
    percent_until_tick_   = curr_tick_size_;
    percent_at_last_tick_ = percent_done;
    update_tick_size();
}

std::string ProgressMeter::position_pad(const GenomicRegion& completed_region) const
{
    assert(position_tab_length_ > 3);
    const auto num_contig_name_letters = completed_region.contig_name().size();
    const auto num_region_end_digits = num_digits(completed_region.end());
    const auto l = num_contig_name_letters + num_region_end_digits + 1; // +1 for ':'
    if (l < position_tab_length_ - 3) {
        return std::string(position_tab_length_ - l - 3, ' ');
    }
    return "";
}

std::string ProgressMeter::completed_pad(const std::string& percent_completed, const std::size_t position_tick_size) const
{
    std::size_t pad {1};
    pad += position_tab_length_ - position_tick_size;
    if (percent_completed.size() < 13) pad += 13 - percent_completed.size();
    return std::string(pad, ' ');
}

std::string ProgressMeter::time_taken_pad(const std::string& time_taken) const
{
    if (time_taken.size() >= 17) return {};
    return std::string(17 - time_taken.size(), ' ');
}

std::string ProgressMeter::ttc_pad(const std::string& ttc) const
{
    if (ttc.size() >= 18) return {};
    return std::string(18 - ttc.size(), ' ');
}

void ProgressMeter::update_tick_size()
{
    constexpr std::size_t max_ticks_to_use {10};
    const auto n = std::min(tick_durations_.size(), max_ticks_to_use);
    if (n > 1) {
        using D = decltype(tick_durations_)::value_type;
        const auto tot = std::accumulate(std::crbegin(tick_durations_), std::next(std::crbegin(tick_durations_), n), D {0});
        using s = std::chrono::seconds;
        const auto mean_tick_duration = std::chrono::duration_cast<s>(tot).count() / n;
        if (mean_tick_duration > 300) {
            curr_tick_size_ = std::max(curr_tick_size_ / 100, min_tick_size_);
        } else if (mean_tick_duration > 60) {
            curr_tick_size_ = std::max(curr_tick_size_ / 10, min_tick_size_);
        } else {
            curr_tick_size_ = max_tick_size_;
        }
    }
}

} // namespace octopus
