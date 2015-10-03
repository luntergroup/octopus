//
//  read_utils.cpp
//  Octopus
//
//  Created by Daniel Cooke on 09/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "read_utils.hpp"

#include <list>
#include <stdexcept>
#include <chrono>
#include <random>

namespace Octopus
{
bool has_minimum_coverage(const std::vector<unsigned>& required_coverage)
{
    return std::all_of(std::cbegin(required_coverage), std::cend(required_coverage),
                       [] (auto coverage) { return coverage == 0; });
}

std::vector<AlignedRead>
sample(const MappableSet<AlignedRead>& reads, const GenomicRegion& region, const unsigned max_coverage,
       const unsigned min_downsample_coverage)
{
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::distance;
    using std::transform;
    
    if (reads.empty()) return std::vector<AlignedRead> {};
    
    auto num_positions = size(region);
    
    auto contained = reads.contained_range(region);
    
    std::vector<unsigned> old_position_coverages = positional_coverage(begin(contained), end(contained), region);
    std::vector<unsigned> required_coverage(num_positions);
    
    transform(cbegin(old_position_coverages), cend(old_position_coverages), begin(required_coverage),
              [min_downsample_coverage] (auto old_coverage) {
                  return std::min(old_coverage, min_downsample_coverage);
              });
    
    std::vector<unsigned> new_position_coverages(num_positions, 0);
    
    auto positions = decompose(region);
    
    static const auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::default_random_engine generator {static_cast<unsigned>(seed)};
    
    // first pass ensures minimum coverage requirements are satisfied
    
    std::vector<AlignedRead> result {};
    result.reserve(reads.count_contained(region) * max_coverage);
    
    std::list<AlignedRead> unsampled_reads(begin(contained), end(contained));
    
    while (!has_minimum_coverage(required_coverage)) {
        std::discrete_distribution<unsigned> covers(cbegin(required_coverage), cend(required_coverage));
        auto sample_position = covers(generator);
        
        auto overlapped = overlap_range(cbegin(unsampled_reads), cend(unsampled_reads), positions[sample_position]);
        
        std::uniform_int_distribution<size_t> read_sampler(0, distance(begin(overlapped), end(overlapped)) - 1);
        
        auto sampled_read_it = std::next(begin(overlapped), read_sampler(generator));
        const AlignedRead& sampled_read {*sampled_read_it};
        result.emplace_back(sampled_read);
        unsampled_reads.erase(sampled_read_it.base());
        
        auto offset    = get_begin(sampled_read) - get_begin(region);
        auto num_reads = size(sampled_read);
        
        transform(begin(required_coverage) + offset, begin(required_coverage) + offset + num_reads,
                  begin(required_coverage) + offset, [] (auto count) { return (count == 0) ? 0 : count - 1; });
        
        transform(begin(new_position_coverages) + offset, begin(new_position_coverages) + offset + num_reads,
                  begin(new_position_coverages) + offset, [] (auto count) { return count + 1; });
    }
    
    //    // second pass increases coverage up to maximum coverage bound
    //    while (!unsampled_reads.empty()) {
    //        // TODO
    //    }
    
    return result;
}

template <typename T>
std::vector<GenomicRegion>
find_target_regions(const T& reads, const GenomicRegion& region, const unsigned max_coverage,
                    const unsigned min_downsample_coverage)
{
    auto above_max_coverage_regions = find_high_coverage_regions(reads, region, max_coverage);
    
    std::vector<GenomicRegion> result {};
    
    if (above_max_coverage_regions.empty()) return result;
    
    result.reserve(above_max_coverage_regions.size());
    
    auto above_min_coverage_regions = find_high_coverage_regions(reads, region, min_downsample_coverage);
    
    std::copy_if(std::cbegin(above_min_coverage_regions), std::cend(above_min_coverage_regions), std::back_inserter(result),
                 [&above_max_coverage_regions] (const auto& r) {
                     return has_contained(std::cbegin(above_max_coverage_regions), std::cend(above_max_coverage_regions), r);
                 });
    
    result.shrink_to_fit();
    
    return result;
}

MappableSet<AlignedRead>
downsample(const MappableSet<AlignedRead>& reads, const unsigned max_coverage, const unsigned min_downsample_coverage)
{
    auto region = get_encompassing_region(std::cbegin(reads), std::cend(reads));
    
    auto regions_to_sample = find_target_regions(reads, region, max_coverage, min_downsample_coverage);
    
    MappableSet<AlignedRead> result {};
    result.reserve(reads.size());
    
    MappableSet<AlignedRead>::const_iterator last_sampled {std::cbegin(reads)};
    
    for (auto& region : regions_to_sample) {
        auto contained = reads.contained_range(region);
        
        if (contained.empty()) continue;
        
        result.insert(last_sampled, std::begin(contained).base());
        
        auto samples = sample(reads, region, max_coverage, min_downsample_coverage);
        
        result.insert(std::make_move_iterator(std::begin(samples)), std::make_move_iterator(std::end(samples)));
        
        last_sampled = std::end(contained).base();
    }
    
    result.insert(last_sampled, std::cend(reads));
    
    result.shrink_to_fit();
    
    return result;
}

// TODO
AlignedRead find_next_segment(const AlignedRead& read, const MappableMap<GenomicRegion::StringType, AlignedRead>& reads)
{
    if (!read.is_chimeric()) {
        throw std::runtime_error {"cannot find next segment as read is not chimeric"};
    }
    
    auto segment_region = read.get_next_segment()->get_inferred_region();
    
    return read;
}

// TODO
MappableSet<AlignedRead> find_chimeras(const AlignedRead& read, const MappableSet<AlignedRead>& reads)
{
    return MappableSet<AlignedRead> {};
}

} // namespace Octopus
