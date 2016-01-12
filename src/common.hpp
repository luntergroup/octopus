//
//  common.hpp
//  Octopus
//
//  Created by Daniel Cooke on 12/05/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_common_hpp
#define Octopus_common_hpp

#include <string>
#include <cstdint>
#include <unordered_map>

#include "genomic_region.hpp"
#include "mappable_set.hpp"
#include "mappable_map.hpp"
#include "aligned_read.hpp"
#include "read_filter.hpp"

namespace Octopus
{
    using ProbabilityType = double;
    using SampleIdType    = std::string;
    
    using SearchRegions = MappableMap<GenomicRegion::ContigNameType, GenomicRegion>;
    
    using ReadContainer = MappableSet<AlignedRead>;
    using ReadMap       = MappableMap<std::string, AlignedRead>;
    
    using ReadFilterer = ReadFilter<ReadContainer::iterator>;
}

#endif
