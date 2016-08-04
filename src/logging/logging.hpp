// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef logging_hpp
#define logging_hpp

#define BOOST_LOG_DYN_LINK 1

#include <iostream>
#include <functional>
#include <sstream>

#include <boost/optional.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/categories.hpp>
#include <boost/iostreams/stream.hpp>
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

namespace octopus { namespace logging {

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
    
    template <typename T> void write(const T& msg) { BOOST_LOG_SEV(lg_, L) << msg; }
    
private:
    src::severity_logger<severity_level> lg_;
};

template <severity_level L, typename T>
Logger<L>& operator<<(Logger<L>& lg, const T& msg)
{
    lg.write(msg);
    return lg;
}

class TraceLogger   : public Logger<severity_level::trace> {};
class DebugLogger   : public Logger<severity_level::debug> {};
class InfoLogger    : public Logger<severity_level::info> {};
class WarningLogger : public Logger<severity_level::warning> {};
class ErrorLogger   : public Logger<severity_level::error> {};
class FatalLogger   : public Logger<severity_level::fatal> {};

//    namespace detail {
//        template <typename CharTp>
//        class BasicLoggerDevice
//        {
//        public:
//            using char_type = CharTp;
//            using category  = boost::iostreams::sink_tag;
//            
//            BasicLoggerDevice() = default;
//            
//            virtual ~BasicLoggerDevice() = default;
//            
//            virtual std::streamsize write(const char_type* s, std::streamsize n) = 0;
//        };
//        
//        template <typename Log>
//        class LoggerBuffer : public BasicLoggerDevice<char>
//        {
//        public:
//            LoggerBuffer() = delete;
//            
//            LoggerBuffer(Log& log, std::streamsize indent_size = 0)
//            : log_ {log}, msg_ {}, newline_indent_ {}
//            {
//                if (indent_size > 0) {
//                    newline_indent_ = '\n' +  std::string(indent_size, ' ');
//                }
//            }
//            
//            LoggerBuffer(const LoggerBuffer&)            = delete;
//            LoggerBuffer& operator=(const LoggerBuffer&) = delete;
//            LoggerBuffer(LoggerBuffer&&)                 = default;
//            LoggerBuffer& operator=(LoggerBuffer&&)      = default;
//            
//            ~LoggerBuffer() noexcept
//            {
//                try {
//                    auto str = msg_.str();
//                    
//                    if (!str.empty() && str.back() == '\n') {
//                        str.pop_back();
//                    }
//                    
//                    if (!newline_indent_.empty()) {
//                        boost::replace_all(str, "\n", newline_indent_);
//                    }
//                    
//                    log_.get() << str;
//                } catch (...) {
//                    return;
//                }
//            }
//            
//            std::streamsize write(const char_type* s, std::streamsize n) override
//            {
//                msg_.write(s, n);
//                return n;
//            }
//            
//        private:
//            std::reference_wrapper<Log> log_;
//            std::ostringstream msg_;
//            std::string newline_indent_;
//        };
//        
//        template <typename Log>
//        class LogStream
//        {
//        public:
//            LogStream(Log& log, std::streamsize indent = 0)
//            : logstream_ {LoggerBuffer<Log> {log, indent}}, os_ {&logstream_}
//            {}
//            
//            LogStream(const LogStream&)            = delete;
//            LogStream& operator=(const LogStream&) = delete;
//            LogStream(LogStream&&)                 = delete;
//            LogStream& operator=(LogStream&&)      = delete;
//            
//            operator std::ostream&() { return os_; }
//            
//        private:
//            boost::iostreams::stream_buffer<LoggerBuffer<Log>> logstream_;
//            std::ostream os_;
//        };
//        
//        template <typename Log>
//        LogStream<Log> make_logstream(Log& log, std::streamsize indent)
//        {
//            return LogStream<Log> {log, indent};
//        }
//    } // namespace detail

template <typename Log>
class LogStream
{
public:
    LogStream() = delete;
    
    LogStream(Log& log, const unsigned indent_size = 0) : log_ {log}, msg_ {}, newline_indent_ {}
    {
        if (indent_size > 0) {
            newline_indent_ = '\n' +  std::string(indent_size, ' ');
        }
    }
    
    LogStream(const LogStream&)            = delete;
    LogStream& operator=(const LogStream&) = delete;
    LogStream(LogStream&&)                 = default;
    LogStream& operator=(LogStream&&)      = default;
    
    ~LogStream()
    {
        auto str = msg_.str();
        
        if (!str.empty() && str.back() == '\n') {
            str.pop_back();
        }
        
        if (!newline_indent_.empty()) {
            boost::replace_all(str, "\n", newline_indent_);
        }
        
        log_.get() << str;
    }
    
    template <typename M> void write(const M& msg) { msg_<< msg; }
    
private:
    std::reference_wrapper<Log> log_;
    std::ostringstream msg_;
    
    std::string newline_indent_;
};

//    template <typename Log>
//    auto logstream(Log& log, std::streamsize indent = 4)
//    {
//        return detail::make_logstream(log, indent);
//    }

template <typename Log>
void log_empty_line(Log& log)
{
    log << "";
}

template <typename Log>
auto stream(Log& log, const unsigned newline_indent = 4)
{
    return LogStream<Log> {log, newline_indent};
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

} // namespace logging
} // namespace octopus

#endif /* logging_hpp */
