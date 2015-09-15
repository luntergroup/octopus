//
//  common.h
//  Octopus
//
//  Created by Daniel Cooke on 12/05/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_common_h
#define Octopus_common_h

#include <string>
#include <cstdint>
#include <unordered_map>

#include "genomic_region.h"
#include "mappable_set.h"
#include "mappable_map.h"
#include "aligned_read.h"

namespace Octopus
{
    using ProbabilityType = double;
    using SampleIdType    = std::string;
    
    using SearchRegions = MappableMap<GenomicRegion::StringType, GenomicRegion>;
    
    using ReadContainer = MappableSet<AlignedRead>;
    using ReadMap       = MappableMap<std::string, AlignedRead>;
}

#endif
