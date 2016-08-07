// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "error_handler.hpp"

#include <string>
#include <vector>
#include <cstddef>

#include <logging/logging.hpp>
#include <utils/string_utils.hpp>

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

void log_error(const Error& error)
{
    logging::ErrorLogger log {};
    
    stream(log) << "A " << error.type() << " error has occured:";
    
    log_empty_line(log);
    
    constexpr unsigned max_line_length {72};
    static const std::string tab {"    "};
    
    auto why_lines = tidy_and_format(error.why(), max_line_length - tab.length());
    
    for (auto& line : why_lines) {
        line.insert(0, tab);
        log << line;
    }
    
    log_empty_line(log);
    
    auto help = std::string {"To help resolve this error "} + error.help();
    
    const auto help_lines = tidy_and_format(help, max_line_length);
    
    for (const auto& line : help_lines) {
        log << line;
    }
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

} // namepace octopus
