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
#include <functional> // std::reference_wrapper

#include "maths.h"

#include <iostream> // TEST

std::vector<unsigned> positional_coverage(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region)
{
    return positional_coverage(reads.cbegin(), reads.cend(), a_region);
}

// Unfortunately the algorithm above is faster than this one. Which is a shame because
// this one is so damn pretty!
//std::vector<unsigned> positional_coverage2(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region)
//{
//    std::vector<unsigned> result(size(a_region), 0);
//    
//    auto positions = decompose(a_region);
//    
//    std::transform(std::cbegin(positions), std::cend(positions), result.begin(),
//                   [&reads] (const auto& position) {
//                       return count_overlapped(reads.cbegin(), reads.cend(), position);
//                   });
//    
//    return result;
//}

unsigned min_coverage(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region)
{
    auto positions_coverage = positional_coverage(reads, a_region);
    return *std::min_element(std::cbegin(positions_coverage), std::cend(positions_coverage));
}

unsigned max_coverage(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region)
{
    auto positions_coverage = positional_coverage(reads, a_region);
    return *std::max_element(std::cbegin(positions_coverage), std::cend(positions_coverage));
}

double mean_coverage(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region)
{
    auto positions_coverage = positional_coverage(reads, a_region);
    return mean(std::cbegin(positions_coverage), std::cend(positions_coverage));
}

double stdev_coverage(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region)
{
    auto positions_coverage = positional_coverage(reads, a_region);
    return stdev(std::cbegin(positions_coverage), std::cend(positions_coverage));
}

std::vector<GenomicRegion> find_high_coverage_regions(const std::vector<AlignedRead>& reads,
                                                      const GenomicRegion& a_region, unsigned maximum_coverage)
{
    std::vector<GenomicRegion> result {};
    
    auto positions_coverage = positional_coverage(reads, a_region);
    
    using Iterator = decltype(positions_coverage)::const_iterator;
    Iterator first {positions_coverage.cbegin()};
    Iterator current {first};
    Iterator last {positions_coverage.cend()};
    Iterator high_range_first, high_range_last;
    GenomicRegion::SizeType high_range_begin, high_range_end;
    
    auto f_is_high_coverage = [maximum_coverage] (unsigned coverage) {
        return coverage > maximum_coverage;
    };
    
    while (current != last) {
        high_range_first = std::find_if(current, last, f_is_high_coverage);
        
        if (high_range_first == last) break;
        
        high_range_last = std::find_if_not(high_range_first, last, f_is_high_coverage);
        
        high_range_begin = get_begin(a_region) + static_cast<GenomicRegion::SizeType>(std::distance(first, high_range_first));
        high_range_end   = high_range_begin + static_cast<GenomicRegion::SizeType>(std::distance(high_range_first, high_range_last));
        
        result.emplace_back(get_contig_name(a_region), high_range_begin, high_range_end);
        
        current = high_range_last;
    }
    
    result.shrink_to_fit();
    
    return result;
}

std::unordered_map<AlignedRead, unsigned> min_coverage_in_read_region(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region)
{
    std::unordered_map<AlignedRead, unsigned> result {};
    result.reserve(reads.size());
    
    auto position_coverages = positional_coverage(reads, a_region);
    auto first_position = get_begin(a_region);
    auto num_positions  = size(a_region);
    
    for (const auto& read : reads) {
        auto first = std::next(std::cbegin(position_coverages), (get_begin(read) <= first_position) ? 0 : get_begin(read) - first_position);
        auto last  = std::next(std::cbegin(position_coverages), std::min(get_end(read) - first_position, num_positions));
        result.emplace(read, *std::min_element(first, last));
    }
    
    return result;
}

std::unordered_map<AlignedRead, unsigned> max_coverage_in_read_region(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region)
{
    std::unordered_map<AlignedRead, unsigned> result {};
    result.reserve(reads.size());
    
    auto position_coverages = positional_coverage(reads, a_region);
    auto first_position = get_begin(a_region);
    auto num_positions  = size(a_region);
    
    for (const auto& read : reads) {
        auto first = std::next(std::cbegin(position_coverages), (get_begin(read) <= first_position) ? 0 : get_begin(read) - first_position);
        auto last  = std::next(std::cbegin(position_coverages), std::min(get_end(read) - first_position, num_positions));
        result.emplace(read, *std::max_element(first, last));
    }
    
    return result;
}

bool has_minimum_coverage(const std::vector<unsigned>& required_coverage)
{
    return std::all_of(required_coverage.cbegin(), required_coverage.cend(),
                       [] (unsigned coverage) {
                           return coverage == 0;
                       });
}

std::vector<GenomicRegion>
find_good_coverage_regions_containing_high_coverage_positions(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region,
                                                              unsigned maximum_coverage, unsigned minimum_downsample_coverage)
{
    auto above_max_coverage_regions = find_high_coverage_regions(reads, a_region, maximum_coverage);
    
    std::vector<GenomicRegion> result {};
    
    if (above_max_coverage_regions.empty()) return result;
    
    result.reserve(above_max_coverage_regions.size());
    
    auto above_min_coverage_regions = find_high_coverage_regions(reads, a_region, minimum_downsample_coverage);
    
    std::copy_if(std::cbegin(above_min_coverage_regions), std::cend(above_min_coverage_regions), std::back_inserter(result),
                 [&above_max_coverage_regions] (const auto& region) {
                     return has_contained(std::cbegin(above_max_coverage_regions), std::cend(above_max_coverage_regions), region);
                 });
    
    result.shrink_to_fit();
    
    return result;
}

std::vector<AlignedRead> sample(std::vector<AlignedRead>::const_iterator first,
                                std::vector<AlignedRead>::const_iterator last,
                                const GenomicRegion& encompassing_region,
                                unsigned maximum_coverage, unsigned minimum_downsample_coverage)
{
    std::cout << "downsampling reads in region " << encompassing_region << std::endl;
    
    if (std::distance(first, last) == 0) return std::vector<AlignedRead> {};
    
    auto num_positions = size(encompassing_region);
    
    std::vector<unsigned> old_position_coverages = positional_coverage(first, last, encompassing_region);
    std::vector<unsigned> required_coverage(num_positions);
    
    std::transform(old_position_coverages.cbegin(), old_position_coverages.cend(), required_coverage.begin(),
                   [minimum_downsample_coverage] (unsigned old_coverage) {
                       return std::min(old_coverage, minimum_downsample_coverage);
                   });
    
    std::vector<unsigned> new_position_coverages(num_positions, 0);
    
    auto positions = decompose(encompassing_region);
    
    static const auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::default_random_engine generator {static_cast<unsigned>(seed)};
    
    // first pass ensures minimum coverage requirements are satisfied
    
    std::vector<AlignedRead> result {};
    result.reserve(std::distance(first, last) * maximum_coverage);
    
    std::list<AlignedRead> unsampled_reads(first, last);
    
    while (!has_minimum_coverage(required_coverage)) {
        std::discrete_distribution<unsigned> covers(required_coverage.cbegin(), required_coverage.cend());
        auto sample_position = covers(generator);
        
        auto overlapped = overlap_range2(unsampled_reads.cbegin(), unsampled_reads.cend(), positions[sample_position]);
        
        std::uniform_int_distribution<std::size_t> read_sampler(0, std::distance(overlapped.begin(), overlapped.end()) - 1);
        
        auto sampled_read_it = std::next(overlapped.begin(), read_sampler(generator));
        const AlignedRead& sampled_read {*sampled_read_it};
        result.emplace_back(sampled_read);
        unsampled_reads.erase(sampled_read_it.base());
        
        auto offset = get_begin(sampled_read) - get_begin(encompassing_region);
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
    
    result.shrink_to_fit();
    
    std::sort(result.begin(), result.end());
    
    return result;
}

std::vector<AlignedRead> downsample(const std::vector<AlignedRead>& reads, unsigned maximum_coverage,
                                    unsigned minimum_downsample_coverage)
{
    auto region = encompassing(reads.cbegin(), reads.cend());
    
    auto regions_to_sample = find_good_coverage_regions_containing_high_coverage_positions(reads, region, maximum_coverage,
                                                                                           minimum_downsample_coverage);
    
    std::vector<AlignedRead> result {};
    result.reserve(reads.size());
    
    std::vector<AlignedRead>::const_iterator last_sampled {reads.cbegin()};
    
    for (auto& region : regions_to_sample) {
        auto contained = contained_range(reads.cbegin(), reads.cend(), region);
        
        if (std::distance(contained.first, contained.second) == 0) continue;
        
        result.insert(result.end(), last_sampled, contained.first);
        
        auto samples = sample(contained.first, contained.second, region, maximum_coverage, minimum_downsample_coverage);
        
        result.insert(result.end(), std::make_move_iterator(std::begin(samples)), std::make_move_iterator(std::end(samples)));
        
        last_sampled = contained.second;
    }
    
    result.insert(result.end(), last_sampled, reads.cend());
    
    result.shrink_to_fit();
    
    return result;
}
