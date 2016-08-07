// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "common.hpp"

namespace octopus {

bool DEBUG_MODE {false}, TRACE_MODE {false};

void log_program_startup()
{
    logging::InfoLogger log {};
    log << "------------------------------------------------------------------------";
    if (TRACE_MODE) {
        stream(log) << "octopus v" << info::Version << " (trace mode)";
    } else if (DEBUG_MODE) {
        stream(log) << "octopus v" << info::Version << " (debug mode)";
    } else {
        stream(log) << "octopus v" << info::Version;
    }
    log << info::CopyrightNotice;
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
