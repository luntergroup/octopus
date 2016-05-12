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

#include "genomic_region.hpp"
#include "read_manager.hpp"
#include "mappable_flat_set.hpp"
#include "mappable_flat_multi_set.hpp"
#include "mappable_map.hpp"
#include "aligned_read.hpp"
#include "read_filter.hpp"

namespace Octopus
{
    extern bool DEBUG_MODE;
    extern bool TRACE_MODE;
    
    const static std::string Octopus_version {"1.0"};
    
    const static std::string Octopus_bug_email {"dcooke@well.ox.ac.uk"};
    
    using SampleIdType = std::string;
    
    using ContigNameType = GenomicRegion::ContigNameType;
    
    using InputRegionMap = MappableSetMap<ContigNameType, GenomicRegion>;
    
    using ReadContainer = MappableFlatMultiSet<AlignedRead>;
    using ReadMap       = MappableMap<SampleIdType, AlignedRead>;
    
    using ReadFilterer = ReadFilter<ReadManager::ReadContainer::iterator>;
}

#endif
