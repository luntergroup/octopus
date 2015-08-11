//
//  vcf_parser.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "vcf_parser.h"

#include <string>

#include "genomic_region.h"
#include "vcf_record.h"

// public methods

VcfParser::VcfParser(const fs::path& file_path)
:
file_path_ {file_path}
{}

VcfHeader VcfParser::fetch_header()
{
    return header_;
}

std::vector<VcfRecord> VcfParser::fetch_records()
{
    std::vector<VcfRecord> result {};
    
    return result;
}

std::vector<VcfRecord> VcfParser::fetch_records(const GenomicRegion& region)
{
    std::vector<VcfRecord> result {};
    
    return result;
}

// private methods

void VcfParser::parse()
{
    
}

// non-member methods

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

//void VcfHeader::insert_line(const std::string& line)
//{
//    if (is_valid_line(line)) {
//        auto content = line.substr(2); // remove ##
//        if (is_format_line(content)) {
//            insert_format_line(content);
//        } else {
//            insert_field_line(content);
//        }
//    } else {
//        throw std::runtime_error {"invalid header line " + line};
//    }
//}

//void insert_format_line(const std::string& line)
//{
//    auto it = std::find(line.cbegin(), line.cend(), '=');
//    formats_.emplace(std::string {line.cbegin(), it}, parse_fields(std::string {std::next(it, 2), std::prev(line.cend())}));
//}
//
//void insert_field_line(const std::string& line)
//{
//    if (!is_valid_field(line)) {
//        throw std::runtime_error {"field line " + line + " has invalid format"};
//    }
//    auto it = std::find(line.cbegin(), line.cend(), '=');
//    fields_.emplace(std::string {std::next(line.cbegin(), 2), it}, std::string {std::next(it), line.cend()});
//}
//
//void insert_lines(const std::string& lines)
//{
//    std::istringstream ss {lines};
//    std::for_each(std::istream_iterator<Line>(ss), std::istream_iterator<Line>(),
//                  [] (const auto& line) { line); });
//}

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
