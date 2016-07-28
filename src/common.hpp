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

#include <boost/optional.hpp>

#include "genomic_region.hpp"
#include "read_manager.hpp"
#include "mappable_flat_set.hpp"
#include "mappable_flat_multi_set.hpp"
#include "mappable_map.hpp"
#include "aligned_read.hpp"
#include "read_filter.hpp"
#include "logging.hpp"

namespace octopus
{
    extern bool DEBUG_MODE;
    extern bool TRACE_MODE;
    
    const static std::string Octopus_version {"1.0"};
    
    const static std::string Octopus_bug_email {"dcooke@well.ox.ac.uk"};
    
    using SampleName = std::string;
    
    using ContigName = GenomicRegion::ContigName;
    
    using InputRegionMap = MappableSetMap<ContigName, GenomicRegion>;
    
    using ReadContainer = MappableFlatMultiSet<AlignedRead>;
    using ReadMap       = MappableMap<SampleName, AlignedRead>;
    
    using ReadFilterer = ReadFilter<ReadManager::ReadContainer::iterator>;
    
    boost::optional<Logging::DebugLogger> get_debug_log();
    boost::optional<Logging::TraceLogger> get_trace_log();
}

#endif
