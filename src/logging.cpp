//
//  logging.cpp
//  Octopus
//
//  Created by Daniel Cooke on 09/03/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "logging.hpp"

#include <iostream>

namespace Octopus
{
namespace Logging
{
    namespace sinks    = boost::log::sinks;
    namespace keywords = boost::log::keywords;
    namespace expr     = boost::log::expressions;
    
    std::ostream& operator<<(std::ostream& os, severity_level level)
    {
        switch (level) {
            case severity_level::trace: os <<  "TRCE"; break;
            case severity_level::debug: os <<  "DEBG"; break;
            case severity_level::info: os <<  "INFO"; break;
            case severity_level::warning: os <<  "WARN"; break;
            case severity_level::error: os <<  "ERRR"; break;
            case severity_level::fatal: os <<  "FATL"; break;
        }
        return os;
    }
    
    void init(boost::optional<boost::filesystem::path> debug_log,
              boost::optional<boost::filesystem::path> trace_log)
    {
        logging::add_console_log
        (
            std::clog,
            keywords::filter =
            (
                severity != severity_level::debug && severity != severity_level::trace
            ),
            keywords::format =
            (
             expr::stream
                << expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "[%Y-%m-%d %H:%M:%S]")
                << " <" << severity
                << "> " << expr::smessage
            )
        );
        
        if (debug_log) {
            logging::add_file_log
            (
             keywords::file_name = debug_log->c_str(),
             keywords::filter =
             (
              severity != severity_level::trace
              ),
             keywords::format =
             (
                expr::stream
                    << expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "[%Y-%m-%d %H:%M:%S]")
                    << " <" << severity
                    << "> " << expr::smessage
              )
             );
        }
        
        if (trace_log) {
            logging::add_file_log
            (
             keywords::file_name = trace_log->c_str(),
             keywords::filter =
             (
              severity != severity_level::debug
              ),
             keywords::format =
             (
              expr::stream
                    << expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "[%Y-%m-%d %H:%M:%S]")
                    << " <" << severity
                    << "> " << expr::smessage
              )
             );
        }
        
        logging::add_common_attributes();
    }

} // namespace Logging
} // namespace Octopus
