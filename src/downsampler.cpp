//
//  downsampler.cpp
//  Octopus
//
//  Created by Daniel Cooke on 19/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "downsampler.hpp"

#include <list>
#include <stdexcept>
#include <random>

#include "read_utils.hpp"

namespace Octopus
{

Downsampler::Downsampler(unsigned max_coverage, unsigned min_coverage)
:
max_coverage_ {max_coverage},
min_coverage_ {min_coverage}
{}

bool has_minimum_coverage(const std::vector<unsigned>& required_coverage)
{
    return std::all_of(std::cbegin(required_coverage), std::cend(required_coverage),
                       [] (auto coverage) { return coverage == 0; });
}

std::vector<AlignedRead>
sample(const MappableSet<AlignedRead>& reads, const GenomicRegion& region,
       const unsigned max_coverage, const unsigned min_coverage)
{
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::distance;
    using std::transform;
    
    if (reads.empty()) return std::vector<AlignedRead> {};
    
    auto num_positions = region_size(region);
    
    auto contained = reads.contained_range(region);
    
    std::vector<unsigned> old_position_coverages {get_positional_coverage(contained, region)};
    std::vector<unsigned> required_coverage(num_positions);
    
    transform(cbegin(old_position_coverages), cend(old_position_coverages), begin(required_coverage),
              [min_coverage] (auto old_coverage) {
                  return std::min(old_coverage, min_coverage);
              });
    
    std::vector<unsigned> new_position_coverages(num_positions, 0);
    
    auto positions = decompose(region);
    
    static std::mt19937_64 generator {1};
    
    // first pass ensures minimum coverage requirements are satisfied
    
    std::vector<AlignedRead> result {};
    result.reserve(reads.count_contained(region) * max_coverage);
    
    std::list<AlignedRead> unsampled_reads {begin(contained), end(contained)};
    
    while (!has_minimum_coverage(required_coverage)) {
        std::discrete_distribution<unsigned> covers {cbegin(required_coverage), cend(required_coverage)};
        
        auto sample_position = covers(generator);
        
        auto overlapped = overlap_range(cbegin(unsampled_reads), cend(unsampled_reads), positions[sample_position]);
        
        std::uniform_int_distribution<size_t> read_sampler(0, distance(begin(overlapped), end(overlapped)) - 1);
        
        auto sampled_read_it = std::next(begin(overlapped), read_sampler(generator));
        const AlignedRead& sampled_read {*sampled_read_it};
        result.emplace_back(sampled_read);
        unsampled_reads.erase(sampled_read_it.base());
        
        const auto offset    = begin_distance(sampled_read, region);
        const auto num_reads = region_size(sampled_read);
        
        transform(begin(required_coverage) + offset, begin(required_coverage) + offset + num_reads,
                  begin(required_coverage) + offset, [] (auto count) { return (count == 0) ? 0 : count - 1; });
        
        transform(begin(new_position_coverages) + offset, begin(new_position_coverages) + offset + num_reads,
                  begin(new_position_coverages) + offset, [] (auto count) { return count + 1; });
    }
    
    // second pass increases coverage up to maximum coverage bound
    
    //    while (!unsampled_reads.empty()) {
    //        // TODO
    //    }
    
    return result;
}

template <typename T>
std::vector<GenomicRegion>
find_target_regions(const T& reads, const GenomicRegion& region, const unsigned max_coverage,
                    const unsigned min_coverage)
{
    auto above_max_coverage_regions = find_high_coverage_regions(reads, region, max_coverage);
    
    std::vector<GenomicRegion> result {};
    
    if (above_max_coverage_regions.empty()) return result;
    
    result.reserve(above_max_coverage_regions.size());
    
    auto above_min_coverage_regions = find_high_coverage_regions(reads, region, min_coverage);
    
    std::copy_if(std::cbegin(above_min_coverage_regions), std::cend(above_min_coverage_regions), std::back_inserter(result),
                 [&above_max_coverage_regions] (const auto& r) {
                     return has_contained(std::cbegin(above_max_coverage_regions), std::cend(above_max_coverage_regions), r);
                 });
    
    result.shrink_to_fit();
    
    return result;
}

MappableSet<AlignedRead>
downsample(const MappableSet<AlignedRead>& reads, const unsigned max_coverage, const unsigned min_coverage)
{
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::make_move_iterator;
    
    if (reads.empty()) return reads;
    
    const auto regions_to_sample = find_target_regions(reads, encompassing_region(reads),
                                                       max_coverage, min_coverage);
    
    if (regions_to_sample.empty()) return reads;
    
    MappableSet<AlignedRead> result {};
    result.reserve(reads.size());
    
    auto last_sampled_itr = cbegin(reads);
    
    for (const auto& region : regions_to_sample) {
        auto contained = reads.contained_range(region);
        
        if (contained.empty()) continue;
        
        result.insert(last_sampled_itr, begin(contained).base());
        
        auto samples = sample(reads, region, max_coverage, min_coverage);
        
        result.insert(make_move_iterator(begin(samples)), make_move_iterator(end(samples)));
        
        last_sampled_itr = end(contained).base();
    }
    
    result.insert(last_sampled_itr, cend(reads));
    
    result.shrink_to_fit();
    
    return result;
}

} // namespace Octopus
