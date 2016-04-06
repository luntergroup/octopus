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

#include <iostream>
#include <functional>
#include <sstream>

#include <boost/optional.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>

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
    
    enum class severity_level { trace, debug, info, warning, error, fatal };
    
    std::ostream& operator<<(std::ostream& os, severity_level level);
    
    BOOST_LOG_INLINE_GLOBAL_LOGGER_DEFAULT(logger, src::severity_logger<severity_level>)
    
    BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity", severity_level)
    
    void init(boost::optional<boost::filesystem::path> debug_log = boost::none,
              boost::optional<boost::filesystem::path> trace_log = boost::none);
    
    template <severity_level L>
    class Logger
    {
    public:
        Logger() : lg_ {logger::get()} {}
        
        template <typename T>
        void write(const T& msg) { BOOST_LOG_SEV(lg_, L) << msg; }
        
    private:
        src::severity_logger<severity_level> lg_;
    };
    
    template <severity_level L, typename T>
    Logger<L>& operator<<(Logger<L>& lg, const T& msg)
    {
        lg.write(msg);
        return lg;
    }
    
    class TraceLogger : public Logger<severity_level::trace> {};
    class DebugLogger : public Logger<severity_level::debug> {};
    class InfoLogger : public Logger<severity_level::info> {};
    class WarningLogger : public Logger<severity_level::warning> {};
    class ErrorLogger : public Logger<severity_level::error> {};
    class FatalLogger : public Logger<severity_level::fatal> {};
    
    template <typename T>
    class LogStream
    {
    public:
        LogStream() = delete;
        
        LogStream(T& log, const unsigned indent_size = 0) : log_ {log}, msg_ {}, indent_ {}
        {
            if (indent_size > 0) {
                indent_ = '\n' +  std::string(indent_size, ' ');
            }
        }
        
        ~LogStream()
        {
            auto str = msg_.str();
            
            if (!str.empty() && str.back() == '\n') {
                str.pop_back();
            }
            
            if (!indent_.empty()) {
                boost::replace_all(str, "\n", indent_);
            }
            
            log_.get() << str;
        }
        
        LogStream(const LogStream&)            = delete;
        LogStream& operator=(const LogStream&) = delete;
        LogStream(LogStream&&)                 = default;
        LogStream& operator=(LogStream&&)      = default;
        
        template <typename M>
        void write(const M& msg) { msg_<< msg; }
        
    private:
        std::reference_wrapper<T> log_;
        std::ostringstream msg_;
        
        std::string indent_;
    };
    
    template <typename T>
    auto stream(T& log, const unsigned indent_size = 4)
    {
        return LogStream<T> {log, indent_size};
    }
    
    inline decltype(auto) stream(std::ostream& os) { return os; }
    
    template <typename T, typename M>
    LogStream<T>& operator<<(LogStream<T>& lg, const M& msg)
    {
        lg.write(msg);
        return lg;
    }
    
    template <typename T, typename M>
    LogStream<T>& operator<<(LogStream<T>&& lg, const M& msg)
    {
        lg.write(msg);
        return lg;
    }
} // namespace Logging
} // namespace Octopus

#endif /* logging_hpp */
