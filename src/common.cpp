//
//  common.cpp
//  Octopus
//
//  Created by Daniel Cooke on 15/03/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "common.hpp"

namespace Octopus
{
    bool DEBUG_MODE {false};
    bool TRACE_MODE {false};
    
    boost::optional<Logging::DebugLogger> get_debug_log()
    {
        static boost::optional<Logging::DebugLogger> log {};
        return DEBUG_MODE ? log : boost::none;
    }
    
    boost::optional<Logging::TraceLogger> get_trace_log()
    {
        static boost::optional<Logging::TraceLogger> log {};
        return DEBUG_MODE ? log : boost::none;
    }
}
