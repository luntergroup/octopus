// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "error_handler.hpp"

#include <string>
#include <vector>
#include <cstddef>
#include <cctype>

#include "exceptions/system_error.hpp"
#include "config/config.hpp"
#include "utils/string_utils.hpp"
#include "logging.hpp"

namespace octopus {

std::string tidy(std::string message)
{
    utils::capitalise_front(message);
    if (!message.empty() && message.back() != '.') message += '.';
    return message;
}

std::vector<std::string> format_as_paragraph(const std::string& message, const std::size_t max_line_length)
{
    if (message.empty()) return {};
    
    auto words = utils::split(message, ' ');
    std::vector<std::string> result {};
    std::string cur_line {};
    
    for (auto& word : words) {
        if (cur_line.length() + word.length() + 1 > max_line_length) {
            result.push_back(std::move(cur_line));
            cur_line = std::move(word);
        } else {
            if (!cur_line.empty()) cur_line += " ";
            cur_line += std::move(word);
        }
    }
    
    if (!cur_line.empty()) result.push_back(std::move(cur_line));
    
    return result;
}

auto tidy_and_format(std::string message, const std::size_t max_line_length)
{
    return format_as_paragraph(tidy(message), max_line_length);
}

void log_error_type(const Error& error, logging::ErrorLogger& log)
{
    auto ls = stream(log);
    const auto type = error.type();
    assert(!type.empty());
    if (type == "unclassified") {
        ls << "An ";
    } else {
        ls << "A ";
    }
    ls << type << " error has occurred:";
}

void log_error_details(const Error& error, logging::ErrorLogger& log)
{
    const auto max_line_length = config::CommandLineWidth;
    static const std::string tab {"    "};
    auto why_lines = tidy_and_format(error.why(), max_line_length - tab.length());
    for (auto& line : why_lines) {
        line.insert(0, tab);
        log << line;
    }
}

void log_error_help(const Error& error, logging::ErrorLogger& log)
{
    std::string help {"To help resolve this error "};
    auto help_detail = error.help();
    assert(!help_detail.empty());
    help_detail.front() = std::tolower(help_detail.front());
    help += help_detail;
    const auto max_line_length = config::CommandLineWidth;
    const auto help_lines = tidy_and_format(help, max_line_length);
    for (const auto& line : help_lines) {
        log << line;
    }
}

void log_error(const Error& error)
{
    logging::ErrorLogger log {};
    log_error_type(error, log);
    log_empty_line(log);
    log_error_details(error, log);
    log_empty_line(log);
    log_error_help(error, log);
}

class BadAlloc : public SystemError
{
    std::string do_where() const override { return "unknown"; }
    std::string do_why() const override { return "system could not satisfy memory request"; }
    std::string do_help() const override
    {
        return "ensure the system sufficient resources or submit an error report";
    }
};

void log_error(const std::bad_alloc& error)
{
    const BadAlloc e {};
    log_error(e);
}

class UnclassifiedError : public Error
{
    std::string do_type() const override { return "unclassified"; }
    std::string do_where() const override { return "unknown"; }
    std::string do_why() const override { return why_; }
    std::string do_help() const override
    {
        return "submit an error report";
    }
    
    std::string why_;
    
public:
    UnclassifiedError(std::string why) : why_ {std::move(why)} {}
};

void log_error(const std::exception& error)
{
    const UnclassifiedError e {error.what()};
    log_error(e);
}

void log_unknown_error()
{
    const UnclassifiedError e {"unknown"};
    log_error(e);
}

} // namespace octopus
