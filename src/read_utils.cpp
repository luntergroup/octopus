//
//  read_utils.cpp
//  Octopus
//
//  Created by Daniel Cooke on 09/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "read_utils.h"

#include <chrono>
#include <random>
#include <list>

std::vector<unsigned> positional_coverage(const MappableSet<AlignedRead>& reads, const GenomicRegion& region)
{
    auto overlapped = reads.overlap_range(region);
    return positional_coverage(overlapped.begin(), overlapped.end(), region);
}

bool has_minimum_coverage(const std::vector<unsigned>& required_coverage)
{
    return std::all_of(required_coverage.cbegin(), required_coverage.cend(),
                       [] (unsigned coverage) {
                           return coverage == 0;
                       });
}

std::vector<AlignedRead> sample(const MappableSet<AlignedRead>& reads, const GenomicRegion& region,
                                unsigned maximum_coverage, unsigned minimum_downsample_coverage)
{
    std::cout << "downsampling reads in region " << region << std::endl;
    
    if (reads.empty()) return std::vector<AlignedRead> {};
    
    auto num_positions = size(region);
    
    auto contained = reads.contained_range(region);
    
    std::vector<unsigned> old_position_coverages = positional_coverage(contained.begin(), contained.end(), region);
    std::vector<unsigned> required_coverage(num_positions);
    
    std::transform(old_position_coverages.cbegin(), old_position_coverages.cend(), required_coverage.begin(),
                   [minimum_downsample_coverage] (unsigned old_coverage) {
                       return std::min(old_coverage, minimum_downsample_coverage);
                   });
    
    std::vector<unsigned> new_position_coverages(num_positions, 0);
    
    auto positions = decompose(region);
    
    static const auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::default_random_engine generator {static_cast<unsigned>(seed)};
    
    // first pass ensures minimum coverage requirements are satisfied
    
    std::vector<AlignedRead> result {};
    result.reserve(reads.count_contained(region) * maximum_coverage);
    
    std::list<AlignedRead> unsampled_reads(contained.begin(), contained.end());
    
    while (!has_minimum_coverage(required_coverage)) {
        std::discrete_distribution<unsigned> covers(required_coverage.cbegin(), required_coverage.cend());
        auto sample_position = covers(generator);
        
        auto overlapped = overlap_range(unsampled_reads.cbegin(), unsampled_reads.cend(), positions[sample_position]);
        
        std::uniform_int_distribution<std::size_t> read_sampler(0, std::distance(overlapped.begin(), overlapped.end()) - 1);
        
        auto sampled_read_it = std::next(overlapped.begin(), read_sampler(generator));
        const AlignedRead& sampled_read {*sampled_read_it};
        result.emplace_back(sampled_read);
        unsampled_reads.erase(sampled_read_it.base());
        
        auto offset = get_begin(sampled_read) - get_begin(region);
        std::transform(required_coverage.begin() + offset,
                       required_coverage.begin() + offset + size(sampled_read),
                       required_coverage.begin() + offset, [] (unsigned count) {
                           return (count == 0) ? 0 : count - 1;
                       });
        std::transform(new_position_coverages.begin() + offset,
                       new_position_coverages.begin() + offset + size(sampled_read),
                       new_position_coverages.begin() + offset, [] (unsigned count) {
                           return count + 1;
                       });
    }
    
    // second pass increases coverage up to maximum coverage bound
    
//    while (!unsampled_reads.empty()) {
//        
//    }
    
    return result;
}

MappableSet<AlignedRead> downsample(const MappableSet<AlignedRead>& reads, unsigned maximum_coverage,
                                    unsigned minimum_downsample_coverage)
{
    auto region = get_encompassing_region(reads.cbegin(), reads.cend());
    
    auto regions_to_sample = find_good_coverage_regions_containing_high_coverage_positions(reads, region, maximum_coverage,
                                                                                           minimum_downsample_coverage);
    
    MappableSet<AlignedRead> result {};
    result.reserve(reads.size());
    
    MappableSet<AlignedRead>::const_iterator last_sampled {reads.cbegin()};
    
    for (auto& region : regions_to_sample) {
        auto contained = reads.contained_range(region);
        
        if (contained.empty()) continue;
        
        result.insert(last_sampled, contained.begin().base());
        
        auto samples = sample(reads, region, maximum_coverage, minimum_downsample_coverage);
        
        result.insert(std::make_move_iterator(std::begin(samples)), std::make_move_iterator(std::end(samples)));
        
        last_sampled = contained.end().base();
    }
    
    result.insert(last_sampled, reads.cend());
    
    result.shrink_to_fit();
    
    return result;
}
