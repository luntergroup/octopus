//
//  read_utils.h
//  Octopus
//
//  Created by Daniel Cooke on 09/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_utils__
#define __Octopus__read_utils__

#include <vector>
#include <algorithm> // std::binary_search

#include "genomic_region.h"
#include "aligned_read.h"

inline bool operator<(const GenomicRegion& lhs, const AlignedRead& rhs)
{
    return lhs < rhs.get_region();
}

inline bool operator<(const AlignedRead& lhs, const GenomicRegion& rhs)
{
    return lhs.get_region() < rhs;
}

template <typename ForwardIterator>
bool has_reads_in_region(ForwardIterator first, ForwardIterator last, const GenomicRegion& a_region)
{
    return std::binary_search(first, last, a_region);
}

inline bool has_reads_in_region(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region)
{
    return has_reads_in_region(reads.cbegin(), reads.cend(), a_region);
}

unsigned get_min_coverage(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region);

unsigned get_max_coverage(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region);

double get_mean_coverage(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region);

#endif /* defined(__Octopus__read_utils__) */
