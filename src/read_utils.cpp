//
//  read_utils.cpp
//  Octopus
//
//  Created by Daniel Cooke on 09/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "read_utils.h"

#include "region_algorithms.h"
#include "maths.h"

std::vector<unsigned> positional_coverage(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region)
{
    auto num_positions = size(a_region);
    
    std::vector<unsigned> result(num_positions, 0);
    
    auto overlapped     = overlap_range(reads.cbegin(), reads.cend(), a_region);
    auto first_position = get_begin(a_region);
    
    std::for_each(overlapped.first, overlapped.second, [&result, first_position, num_positions] (const auto& read) {
        auto first = std::next(result.begin(), (get_begin(read) <= first_position) ? 0 : get_begin(read) - first_position);
        auto last  = std::next(result.begin(), std::min(get_end(read) - first_position, num_positions));
        std::transform(first, last, first, [] (unsigned count) { return count + 1; });
    });
    
    return result;
}

// Unfortunatly the algorithm above is faster than this one. Which is a shame because
// this one is so damn pretty.
//std::vector<unsigned> positional_coverage(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region)
//{
//    auto num_positions = size(a_region);
//    
//    std::vector<unsigned> result(num_positions, 0);
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
                                                      const GenomicRegion& a_region, unsigned max_positional_coverage)
{
    auto positions_coverage = positional_coverage(reads, a_region);
    
    std::vector<GenomicRegion> result {};
    
    auto first   = positions_coverage.cbegin();
    auto current = first;
    auto last    = positions_coverage.cend();
    
    decltype(positions_coverage)::const_iterator first_high, last_high;
    GenomicRegion::SizeType high_begin, high_end;
    
    while (current != last) {
        first_high = std::find_if(current, last, [max_positional_coverage] (unsigned coverage) {
            return coverage > max_positional_coverage;
        });
        
        if (first_high == last) break;
        
        last_high = std::find_if_not(first_high, last, [max_positional_coverage] (unsigned coverage) {
            return coverage > max_positional_coverage;
        });
        
        high_begin = get_begin(a_region) + static_cast<GenomicRegion::SizeType>(std::distance(first, first_high));
        high_end   = get_begin(a_region) + static_cast<GenomicRegion::SizeType>(std::distance(first, last_high));
        
        result.emplace_back(get_contig_name(a_region), high_begin, high_end);
        
        current = last_high;
    }
    
    result.shrink_to_fit();
    
    return result;
}

std::vector<AlignedRead> downsample(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region,
                                    unsigned max_positional_coverage)
{
    auto region = get_encompassing(reads.front(), *rightmost_mappable(reads.cbegin(), reads.cend()));
    
    if (max_coverage(reads, region) <= max_positional_coverage) return reads;
    
    return reads;
}
