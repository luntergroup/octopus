// Copyright (c) 2015-2020 Daniel Cooke
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

} // namespace detail

void log_program_startup()
{
    logging::InfoLogger log {};
    const auto banner = detail::make_banner();
    log << banner;
    std::ostringstream ss {};
    ss << "octopus v" << config::Version;
    log << ss.str();
    log << config::CopyrightNotice;
    log << banner;
}

void log_command_line_options(const options::OptionMap& options)
{
    static auto debug_log = logging::get_debug_log();
    if (debug_log) {
        *debug_log << "Program options: ";
        *debug_log << options::to_string(options);
    }
}

} // namespace octopus
