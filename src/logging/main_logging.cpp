// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "main_logging.hpp"

#include <sstream>

#include "config/config.hpp"
#include "config/common.hpp"

namespace octopus {

namespace detail {
    std::string make_banner()
    {
        return std::string(config::CommandLineWidth, '-');
    }
}

void log_program_startup()
{
    logging::InfoLogger log {};
    
    const auto banner = detail::make_banner();
    
    log << banner;
    
    std::ostringstream ss {};
    
    ss << "octopus v" << config::Version;
    
    if (TRACE_MODE) {
        ss << " (trace mode)";
    } else if (DEBUG_MODE) {
        ss << " (debug mode)";
    }
    
    log << ss.str();
    
    log << config::CopyrightNotice;
    log << banner;
}

} // namespace octopus
