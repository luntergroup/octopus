//
//  common.cpp
//  octopus
//
//  Created by Daniel Cooke on 15/03/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "common.hpp"

namespace octopus {

bool DEBUG_MODE {false}, TRACE_MODE {false};

void log_program_startup()
{
    logging::InfoLogger log {};
    log << "------------------------------------------------------------------------";
    if (TRACE_MODE) {
        stream(log) << "octopus v" << info::VERSION << " (trace mode)";
    } else if (DEBUG_MODE) {
        stream(log) << "octopus v" << info::VERSION << " (debug mode)";
    } else {
        stream(log) << "octopus v" << info::VERSION;
    }
    log << info::COPYRIGHT_NOTICE;
    log << "------------------------------------------------------------------------";
}

namespace logging {

boost::optional<DebugLogger> get_debug_log()
{
    if (DEBUG_MODE) {
        return DebugLogger {};
    }
    return boost::none;
}

boost::optional<TraceLogger> get_trace_log()
{
    if (TRACE_MODE) {
        return TraceLogger {};
    }
    return boost::none;
}

} // namespace debug
    
} // namespace octopus
