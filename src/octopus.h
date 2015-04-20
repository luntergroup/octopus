//
//  octopus.h
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__octopus__
#define __Octopus__octopus__

#include <vector>

#include "read_manager.h"

class GenomicRegion;
class AlignedRead;
class Variant;

using Variants = std::vector<Variant>;

void run_octopus();

GenomicRegion next_sub_region(const GenomicRegion& the_search_region, const GenomicRegion& previous_sub_region,
                              const ReadManager::SampleReadMap& the_reads, const Variants& the_candidates,
                              unsigned max_variants_in_region, unsigned max_region_size, unsigned max_region_overlap);

#endif /* defined(__Octopus__octopus__) */
