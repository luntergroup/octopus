//
//  logging.hpp
//  Octopus
//
//  Created by Daniel Cooke on 09/03/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef logging_hpp
#define logging_hpp

#define BOOST_LOG_DYN_LINK 1

#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/sources/global_logger_storage.hpp>
#include <boost/move/utility.hpp>
#include <boost/log/sources/logger.hpp>

namespace Octopus
{
namespace Logging
{
    namespace logging  = boost::log;
    namespace src      = boost::log::sources;
    
    BOOST_LOG_INLINE_GLOBAL_LOGGER_DEFAULT(logger, src::severity_logger<logging::trivial::severity_level>)
    
    void init(boost::optional<boost::filesystem::path> log);
} // namespace Logging
} // namespace Octopus

#endif /* logging_hpp */
