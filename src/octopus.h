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

GenomicRegion next_sub_region(const GenomicRegion& the_search_region, const GenomicRegion& the_previous_sub_region,
                              const ReadManager::SampleReadMap& the_reads, const Variants& the_candidates,
                              unsigned max_variants, unsigned max_indicators);

#endif /* defined(__Octopus__octopus__) */
