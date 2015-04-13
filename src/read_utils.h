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

#include "aligned_read.h"

template <typename ForwardIterator>
unsigned get_min_coverage(ForwardIterator first, ForwardIterator last)
{
    return 0;
}

unsigned get_min_coverage(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region);

unsigned get_max_coverage(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region);

double get_mean_coverage(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region);

#endif /* defined(__Octopus__read_utils__) */
