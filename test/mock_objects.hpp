//
//  mock_objects.h
//  Octopus
//
//  Created by Daniel Cooke on 22/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_mock_objects_h
#define Octopus_mock_objects_h

#include <vector>
#include <algorithm>
#include <chrono>
#include <random>
#include <iterator>
#include <cstddef>

#include <basics/genomic_region.hpp>
#include <basics/aligned_read.hpp>
#include <basics/cigar_string.hpp>

namespace octopus { namespace test {

inline auto generate_random_regions(const GenomicRegion::Position contig_size,
                                    const GenomicRegion::Position mean_region_size,
                                    const std::size_t num_regions)
{
    std::vector<GenomicRegion> result {};
    result.reserve(num_regions);
    
    GenomicRegion::ContigName contig_name = "test";
    
    const auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    
    std::default_random_engine generator {static_cast<unsigned>(seed)};
    
    std::uniform_int_distribution<GenomicRegion::Position> begin_dist(0, contig_size - mean_region_size - 1);
    std::geometric_distribution<GenomicRegion::Position> size_dist(1.0 / mean_region_size);
    
    std::generate_n(std::back_inserter(result), num_regions,
                    [&] () {
                        auto begin = begin_dist(generator);
                        auto size  = size_dist(generator);
                        if (begin + size >= contig_size) begin = contig_size - size;
                        return GenomicRegion {contig_name, begin, begin + size};
                    });
    
    std::sort(std::begin(result), std::end(result));
    
    return result;
}

inline auto get_mock_region()
{
    return GenomicRegion {"1", 0, 1};
}

inline auto region_to_read(const GenomicRegion& region)
{
    return AlignedRead {};
}

template <typename Container>
auto regions_to_reads(const Container& regions)
{
    std::vector<AlignedRead> result {};
    result.reserve(regions.size());
    for (const auto& region : regions) {
        result.push_back(region_to_read(region));
    }
    return result;
}

inline auto get_mock_aligned_read(std::string sequence)
{
    return AlignedRead {};
}

} // namespace test
} // namespace octopus

#endif
