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
            case severity_level::debug: os <<  "DEBG"; break;
            case severity_level::info: os <<  "INFO"; break;
            case severity_level::warning: os <<  "WARN"; break;
            case severity_level::error: os <<  "ERRR"; break;
            case severity_level::fatal: os <<  "FATL"; break;
        }
        return os;
    }
    
    void init(boost::optional<boost::filesystem::path> log)
    {
        logging::add_console_log
        (
            std::clog,
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
        
        if (log) {
            logging::add_file_log
            (
             keywords::file_name = log->c_str(),
             keywords::rotation_size = 10 * 1024 * 1024,
             keywords::time_based_rotation = sinks::file::rotation_at_time_point(0, 0, 0),
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
