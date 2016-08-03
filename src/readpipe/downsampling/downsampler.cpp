//
//  downsampler.cpp
//  Octopus
//
//  Created by Daniel Cooke on 19/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "downsampler.hpp"

#include <vector>
#include <deque>
#include <functional>
#include <iterator>
#include <algorithm>
#include <random>
#include <cassert>

#include <concepts/mappable.hpp>
#include <basics/mappable_reference_wrapper.hpp>
#include <concepts/mappable_range.hpp>
#include <utils/mappable_algorithms.hpp>
#include <utils/read_algorithms.hpp>

#include "timers.hpp"

namespace octopus { namespace readpipe {

namespace {

using ReadWrapper = MappableReferenceWrapper<const AlignedRead>;

using PositionCoverages = std::vector<unsigned>;

bool has_minimum_coverage(const PositionCoverages& required_coverage)
{
    return std::all_of(std::cbegin(required_coverage), std::cend(required_coverage),
                       [] (const auto coverage) { return coverage == 0; });
}

template <typename ForwardIt>
auto calculate_minimum_coverages(const ForwardIt first, const ForwardIt last,
                                 const GenomicRegion& region, const unsigned min_coverage)
{
    auto result = calculate_positional_coverage(first, last, region);
    
    std::transform(std::cbegin(result), std::cend(result), std::begin(result),
                   [min_coverage] (const auto coverage) {
                       return std::min(coverage, min_coverage);
                   });
    
    return result;
}

auto sample(const PositionCoverages& required_coverage)
{
    static std::default_random_engine generator {};
    
    // TODO: Do we really need to keep regenerating this distribution?
    std::discrete_distribution<std::size_t> dist {
        std::cbegin(required_coverage), std::cend(required_coverage)
    };
    
    return dist(generator);
}

template <typename BidirIt>
BidirIt random_sample(const BidirIt first, const BidirIt last)
{
    static std::default_random_engine generator {};
    
    std::uniform_int_distribution<std::size_t> dist(0, std::distance(first, last) - 1);
    
    return std::next(first, dist(generator));
}

template <typename T>
auto random_sample(const OverlapRange<T>& range)
{
    return random_sample(std::begin(range), std::end(range)).base();
}

template <typename BidirIt>
auto pick_sample(BidirIt first_unsampled, BidirIt last_unsampled,
                 const std::vector<GenomicRegion>& positions,
                 const PositionCoverages& required_coverage,
                 const AlignedRead::RegionType::Size max_read_size)
{
    assert(first_unsampled < last_unsampled);
    
    const auto candidates = overlap_range(first_unsampled, last_unsampled,
                                          positions[sample(required_coverage)],
                                          max_read_size);
    
    assert(!candidates.empty());
    
    return random_sample(candidates);
}

void reduce(PositionCoverages& coverages, const ReadWrapper& read, const GenomicRegion& region)
{
    assert(!begins_before(read, region));
    
    const auto read_offset = begin_distance(region, read);
    
    const auto first = std::next(std::begin(coverages), read_offset);
    
    auto last = std::next(first, region_size(read));
    
    if (last > std::end(coverages)) last = std::end(coverages);
    
    std::transform(first, last, first, [] (const auto x) { return (x > 0) ? x - 1 : 0; });
}

template <typename ForwardIt>
ForwardIt shift_left(ForwardIt first, ForwardIt last)
{
    return std::rotate(first, std::next(first), last);
}

template <typename BidirIt>
BidirIt shift_right(BidirIt first, BidirIt last)
{
    const auto rfirst = std::make_reverse_iterator(last);
    const auto rlast  = std::make_reverse_iterator(first);
    return std::rotate(rfirst, std::next(rfirst), rlast).base();
}

template <typename BidirIt>
void remove_sample(BidirIt& first, BidirIt sample, BidirIt& last)
{
    assert(first <= sample && sample < last);
    // rotating is linear time so it's cheaper to shift the sample into the closer of the
    // two sampled groups.
    if (std::distance(first, sample) < std::distance(sample, last)) {
        if (sample != first) {
            shift_right(first, std::next(sample));
        }
        ++first;
    } else {
        if (sample != std::prev(last)) {
            shift_left(sample, last);
        }
        --last;
    }
}

auto extract_sampled(std::vector<ReadWrapper>& reads,
                     std::vector<ReadWrapper>::iterator first_unsampled,
                     std::vector<ReadWrapper>::iterator last_unsampled)
{
    reads.erase(first_unsampled, last_unsampled);
    
    std::sort(std::begin(reads), std::end(reads));
    
    std::vector<AlignedRead> result {};
    result.reserve(reads.size());
    
    std::transform(std::begin(reads), std::end(reads), std::back_inserter(result),
                   [] (const auto& wrapper) -> AlignedRead { return wrapper.get(); });
    
    return result;
}

} // namespace

template <typename InputIt>
auto sample(const InputIt first_read, const InputIt last_read, const GenomicRegion& region,
            const unsigned target_coverage)
{
    if (first_read == last_read) return std::vector<AlignedRead> {};
    
    const auto positions = decompose(region);
    
    auto required_coverage = calculate_minimum_coverages(first_read, last_read, region,
                                                         target_coverage);
    
    assert(positions.size() == required_coverage.size());
    
    std::vector<ReadWrapper> reads {first_read, last_read};
    
    const auto max_read_size = size(largest_region(reads)); // for efficient overlap detection
    
    // The reads are partitioned into three groups: sampled | unsampled | sampled
    // which allows a sampled read to be optimally moved out of the unsampled partition
    auto first_unsampled_itr = std::begin(reads);
    auto last_unsampled_itr  = std::end(reads);
    
    while (!has_minimum_coverage(required_coverage)) {
        const auto sampled_itr = pick_sample(first_unsampled_itr, last_unsampled_itr,
                                             positions, required_coverage, max_read_size);
        
        reduce(required_coverage, *sampled_itr, region);
        
        remove_sample(first_unsampled_itr, sampled_itr, last_unsampled_itr);
    }
    
    return extract_sampled(reads, first_unsampled_itr, last_unsampled_itr);
}

namespace {

// Look for regions with coverage above max_coverage interconnected by positions
// with coverage abouve min_coverage. The idea being samples are better taken from the
// larger joined regions.
auto find_target_regions(const ReadContainer& reads, const unsigned max_coverage, const unsigned min_coverage)
{
    const auto above_max_coverage_regions = find_high_coverage_regions(reads, max_coverage);
    
    std::vector<GenomicRegion> result {};
    
    if (above_max_coverage_regions.empty()) return result;
    
    result.reserve(above_max_coverage_regions.size());
    
    auto above_min_coverage_regions = find_high_coverage_regions(reads, min_coverage);
    
    std::copy_if(std::make_move_iterator(std::begin(above_min_coverage_regions)),
                 std::make_move_iterator(std::end(above_min_coverage_regions)),
                 std::back_inserter(result),
                 [&above_max_coverage_regions] (const auto& region) {
                     return has_contained(std::cbegin(above_max_coverage_regions),
                                          std::cend(above_max_coverage_regions),
                                          region);
                 });
    
    return result;
}

void prepend(std::vector<AlignedRead>&& src, std::deque<AlignedRead>& dest)
{
    dest.insert(std::begin(dest),
                std::make_move_iterator(std::begin(src)),
                std::make_move_iterator(std::end(src)));
}

} // namespace

std::size_t sample(ReadContainer& reads, const unsigned trigger_coverage, const unsigned target_coverage)
{
    using std::begin; using std::end; using std::make_move_iterator;
    
    const auto prior_num_reads = reads.size();
    
    if (prior_num_reads == 0) return 0;
    
    const auto targets = find_target_regions(reads, trigger_coverage, target_coverage);
    
    if (targets.empty()) return 0;
    
    std::deque<AlignedRead> sampled {};
    
    // Downsample in reverse order because erasing near back of MappableFlatMultiSet is much
    // cheaper than erasing near front.
    std::for_each(std::crbegin(targets), std::crend(targets),
                  [&] (const auto& region) {
                      const auto contained = bases(contained_range(begin(reads), end(reads), region));
                      
                      prepend(sample(begin(contained), end(contained), region, target_coverage), sampled);
                      
                      reads.erase(begin(contained), end(contained));
                  });
    
    reads.insert(make_move_iterator(begin(sampled)), make_move_iterator(end(sampled)));
    
    return prior_num_reads - reads.size();
}

// Downsampler

Downsampler::Downsampler(unsigned trigger_coverage, unsigned target_coverage)
:
trigger_coverage_ {trigger_coverage},
target_coverage_ {target_coverage}
{
    if (target_coverage > trigger_coverage) {
        target_coverage_ = trigger_coverage;
    }
}

std::size_t Downsampler::downsample(ReadContainer& reads) const
{
    return sample(reads, trigger_coverage_, target_coverage_);
}

} // namespace readpipe
} // namespace octopus
