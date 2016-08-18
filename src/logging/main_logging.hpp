// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef main_logging_hpp
#define main_logging_hpp

#include <string>

#include "logging.hpp"

namespace octopus {

namespace detail {
    std::string make_banner();
}

void log_program_startup(); // Always uses InfoLogger

template <typename Log>
void log_program_end(Log& log)
{
    log << detail::make_banner();
}

inline void log_program_end()
{
    logging::InfoLogger log {};
    log_program_end(log);
}
    
} // namespace octopus

#endif
