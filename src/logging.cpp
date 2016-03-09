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
    
    void init(boost::optional<boost::filesystem::path> log)
    {
        logging::add_console_log
        (
            std::cout,
            keywords::filter =
            (
                logging::trivial::severity >= logging::trivial::info
            ),
            keywords::format =
            (
             expr::stream
                << expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "[%Y-%m-%d %H:%M:%S]")
                << " <" << logging::trivial::severity
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
                    << " <" << logging::trivial::severity
                    << "> " << expr::smessage
              )
             );
        }
        
        logging::add_common_attributes();
    }

} // namespace Logging
} // namespace Octopus
