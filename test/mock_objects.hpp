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

inline std::vector<GenomicRegion> generate_random_regions(GenomicRegion::SizeType contig_size,
                                                          GenomicRegion::SizeType mean_region_size,
                                                          std::size_t num_regions)
{
    std::vector<GenomicRegion> result {};
    result.reserve(num_regions);
    
    GenomicRegion::ContigNameType contig_name = "test";
    
    const auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator {static_cast<unsigned>(seed)};
    
    std::uniform_int_distribution<GenomicRegion::SizeType> begin_dist(0, contig_size - mean_region_size - 1);
    std::geometric_distribution<GenomicRegion::SizeType> size_dist(1.0 / mean_region_size);
    
    std::generate_n(std::back_inserter(result), num_regions,
                    [&contig_name, contig_size, mean_region_size, &generator, &begin_dist, &size_dist] () {
                        auto begin = begin_dist(generator);
                        auto size  = size_dist(generator);
                        if (begin + size >= contig_size) begin = contig_size - size;
                        return GenomicRegion {contig_name, begin, begin + size};
                    });
    
    std::sort(result.begin(), result.end());
    
    return result;
}


inline GenomicRegion get_mock_region()
{
    return GenomicRegion {"1", 0, 1};
}

inline AlignedRead get_mock_aligned_read(std::string the_sequence)
{
    return AlignedRead {
        get_mock_region(),
        std::move(the_sequence),
        AlignedRead::Qualities {},
        CigarString {},
        0,
        AlignedRead::Flags {}
    };
}

#endif
