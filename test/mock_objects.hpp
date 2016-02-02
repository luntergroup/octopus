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

#include "genomic_region.hpp"
#include "aligned_read.hpp"
#include "cigar_string.hpp"
#include "aligned_read_builder.hpp"

inline auto generate_random_regions(const GenomicRegion::SizeType contig_size,
                                    const GenomicRegion::SizeType mean_region_size,
                                    const std::size_t num_regions)
{
    std::vector<GenomicRegion> result {};
    result.reserve(num_regions);
    
    GenomicRegion::ContigNameType contig_name = "test";
    
    const auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    
    std::default_random_engine generator {static_cast<unsigned>(seed)};
    
    std::uniform_int_distribution<GenomicRegion::SizeType> begin_dist(0, contig_size - mean_region_size - 1);
    std::geometric_distribution<GenomicRegion::SizeType> size_dist(1.0 / mean_region_size);
    
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
    return AlignedReadBuilder().set_region(region).set_sequence("").build_once();
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
    return AlignedReadBuilder().set_region(get_mock_region()).set_sequence(std::move(sequence)).build_once();
}

#endif
