//
//  vcf_header.cpp
//  Octopus
//
//  Created by Daniel Cooke on 28/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "vcf_header.h"

#include <istream>
#include <sstream>
#include <algorithm> // std::for_each, std::count, std::find, std::transform
#include <utility>   // std::move
#include <iterator>  // std::cbegin, std::cend, std::next, std::prev

#include <iostream> // TEST

bool is_valid_line(const std::string& line);
bool is_valid_field(const std::string& str);
bool is_format_line(const std::string& line);
std::unordered_map<std::string, std::string> parse_fields(const std::string& fields);

class UnknownKey : std::runtime_error {
public:
    UnknownKey(std::string key)
    :
    runtime_error {"key is not in header"},
    key_ {std::move(key)}
    {}
    
    const char* what() const noexcept
    {
        return (std::string{runtime_error::what()} + ": " + key_).c_str();
    }
    
private:
    std::string key_;
};

// public methods

VcfHeader::VcfHeader(const std::string& lines)
{
    insert_lines(lines);
}

void VcfHeader::insert_line(const std::string& line)
{
    if (is_valid_line(line)) {
        auto content = line.substr(2); // remove ##
        if (is_format_line(content)) {
            insert_format_line(content);
        } else {
            insert_field_line(content);
        }
    } else {
        throw std::runtime_error {"invalid header line " + line};
    }
}

const std::string& VcfHeader::get_file_format() const noexcept
{
    return file_format_;
}

unsigned VcfHeader::get_num_samples() const noexcept
{
    return static_cast<unsigned>(samples_.size());
}

std::vector<std::string> VcfHeader::get_samples() const
{
    return samples_;
}

unsigned VcfHeader::get_field_cardinality(const std::string& key, const VcfRecord& record) const
{
    return 0;
}

VcfType VcfHeader::get_typed_value(const std::string& format_key, const std::string& field_key,
                                   const std::string& value) const
{
    if (formats_.count(format_key) == 0) throw UnknownKey {format_key};
    auto er = formats_.equal_range(format_key);
    auto has_key = [&field_key] (const auto& p) { return p.second.at("ID") == field_key; };
    auto it = std::find_if(er.first, er.second, has_key);
    if (it == er.second) throw UnknownKey {field_key};
    return vcf_type_factory(it->second.at("Type"), value);
}

VcfType VcfHeader::get_typed_value(const std::string& field_key, const std::string& value) const
{
    return vcf_type_factory("String", value);
}

// private methods

template <char Delim>
struct Token
{
    std::string data;
    
    operator std::string() const
    {
        return data;
    }
};

template <char Delim>
std::istream& operator>>(std::istream& str, Token<Delim>& data)
{
    std::getline(str, data.data, Delim);
    return str;
}

using Line  = Token<'\n'>;
using Field = Token<','>;

void VcfHeader::insert_format_line(const std::string& line)
{
    auto it = std::find(line.cbegin(), line.cend(), '=');
    formats_.emplace(std::string {line.cbegin(), it}, parse_fields(std::string {std::next(it, 2), std::prev(line.cend())}));
}

void VcfHeader::insert_field_line(const std::string& line)
{
    if (!is_valid_field(line)) {
        throw std::runtime_error {"field line " + line + " has invalid format"};
    }
    auto it = std::find(line.cbegin(), line.cend(), '=');
    fields_.emplace(std::string {std::next(line.cbegin(), 2), it}, std::string {std::next(it), line.cend()});
}

void VcfHeader::insert_lines(const std::string& lines)
{
    std::istringstream ss {lines};
    std::for_each(std::istream_iterator<Line>(ss), std::istream_iterator<Line>(),
                  [this] (const auto& line) { this->insert_line(line); });
}

// non-member methods

bool is_valid_line(const std::string& line)
{
    return line.size() > 3 && line[0] == '#' && line[1] == '#' && std::find(line.cbegin(), line.cend(), '=') != line.cend();
}

bool is_valid_field(const std::string& str)
{
    return std::count(std::cbegin(str), std::cend(str), '=') == 1; // anything else?
}

bool is_format_line(const std::string& line)
{
    return *std::next(std::find(line.cbegin(), line.cend(), '=')) == '<' && line.back() == '>';
}

std::pair<std::string, std::string> parse_field(const std::string& field)
{
    auto it = std::find(field.cbegin(), field.cend(), '=');
    return std::make_pair(std::string {field.cbegin(), it}, std::string {std::next(it), field.cend()});
}

std::unordered_map<std::string, std::string> parse_fields(const std::string& fields)
{
    std::istringstream ss {fields};
    std::unordered_map<std::string, std::string> result {};
    std::transform(std::istream_iterator<Field>(ss), std::istream_iterator<Field>(), std::inserter(result, result.begin()),
                  [] (const auto& field) { return parse_field(field); });
    return result;
}

std::ostream& operator<<(std::ostream& os, const std::unordered_map<std::string, std::string>& fields)
{
    os << "<";
    for (const auto& p : fields) {
        os << p.first << "=" << p.second << ",";
    }
    os << ">";
    return os;
}

std::ostream& operator<<(std::ostream& os, const VcfHeader& header)
{
    os << "##" << header.file_format_ << std::endl;
    
    for (const auto& field : header.fields_) {
        os << "##" << field.first << "=" << field.second << std::endl;
    }
    
    for (const auto& format : header.formats_) {
        os << "##" << format.first << "=" << format.second << std::endl;
    }
    
    return os;
}
