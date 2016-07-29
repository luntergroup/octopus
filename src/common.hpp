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
#include <vector>
#include <cstdint>

#include <boost/optional.hpp>

#include "genomic_region.hpp"
#include "read_manager.hpp"
#include "mappable_flat_set.hpp"
#include "mappable_flat_multi_set.hpp"
#include "mappable_map.hpp"
#include "aligned_read.hpp"
#include "logging.hpp"

namespace octopus
{
    extern bool DEBUG_MODE;
    extern bool TRACE_MODE;
    
    namespace info
    {
        const static std::string VERSION {"1.0"};
        const static std::string BUG_EMAIL {"dcooke@well.ox.ac.uk"};
        const static std::vector<std::string> AUTHORS {"Daniel Cooke"};
        const static std::string COPYRIGHT_NOTICE {"Copyright (c) 2016 University of Oxford"};
    }
    
    using SampleName = std::string;
    
    using ContigName = GenomicRegion::ContigName;
    
    using InputRegionMap = MappableSetMap<ContigName, GenomicRegion>;
    
    using ReadContainer = MappableFlatMultiSet<AlignedRead>;
    using ReadMap       = MappableMap<SampleName, AlignedRead>;
    
    namespace logging
    {
        boost::optional<DebugLogger> get_debug_log();
        boost::optional<TraceLogger> get_trace_log();
    }
}

#endif
