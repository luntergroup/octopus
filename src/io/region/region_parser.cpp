// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "region_parser.hpp"

#include <regex>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <functional>

#include <exceptions/user_error.hpp>
#include <utils/string_utils.hpp>

namespace octopus { namespace io {

class MalformedRegion : public UserError
{
    std::string do_where() const override { return "parse_region"; }
    
    std::string do_why() const override
    {
        std::ostringstream ss {};
        ss << "the region you specified '";
        ss << region_;
        ss << "' ";
        ss << why_;
        return ss.str();
    }
    
    std::string do_help() const override
    {
        return "ensure the region is in the format contig[:begin][-][end]";
    }
    
    std::string region_;
    
public:
    MalformedRegion(std::string region) : region_ {std::move(region)} {}
    
    MalformedRegion(std::string region, std::string why)
    : region_ {std::move(region)}, why_ {std::move(why)} {}
    
    std::string why_ = "is not formatted correctly";
};

class MissingReferenceContig : public UserError
{
    std::string do_where() const override { return "parse_region"; }
    
    std::string do_why() const override
    {
        return "the region you specified '" + region_ + "' refers to a contig not in the reference " + reference_;
    }
    
    std::string do_help() const override
    {
        return "ensure the region is in the format contig[:begin][-][end], where contig refers to"
        " a contig in the reference " + reference_;
    }
    
    std::string region_, reference_;
    
public:
    MissingReferenceContig(std::string region, std::string reference)
    : region_ {std::move(region)}, reference_ {std::move(reference)} {}
};

GenomicRegion parse_region(std::string region, const ReferenceGenome& reference)
{
    using Position = GenomicRegion::Position;
    
    region.erase(std::remove(std::begin(region), std::end(region), ','), std::end(region));
    
    static const std::regex re {"([^:]+)(?::(\\d+)(-)?(\\d*))?"};
    
    std::smatch match;
    
    if (std::regex_match(region, match, re) && match.size() == 5) {
        GenomicRegion::ContigName contig {match.str(1)};
        
        if (!reference.has_contig(contig)) {
            throw MissingReferenceContig {region, reference.name()};
        }
        
        const auto contig_size = reference.contig_size(contig);
        
        Position begin {0}, end {0};
        
        if (match.length(2) == 0) {
            end = contig_size;
        } else {
            begin = static_cast<Position>(std::stoul(match.str(2)));
            
            if (match.length(3) == 0) {
                end = begin + 1;
            } else if (match.str(4).empty()) {
                end = contig_size;
            } else {
                end = static_cast<Position>(std::stoul(match.str(4)));
            }
            
            if (begin > end) {
                throw MalformedRegion {region, "has begin greater than end"};
            }
            
            if (begin > contig_size) {
                begin = contig_size;
            }
            
            if (end > contig_size) end = contig_size;
        }
        
        return GenomicRegion {std::move(contig), begin, end};
    }
    
    throw MalformedRegion {region};
}

namespace {

struct Line
{
    std::string line_data;
    
    operator std::string() const
    {
        return line_data;
    }
};

std::istream& operator>>(std::istream& is, Line& data)
{
    std::getline(is, data.line_data);
    
    if (!data.line_data.empty() && data.line_data.back() == '\r') {
        data.line_data.pop_back();
    }
    
    return is;
}

bool is_bed_file(const boost::filesystem::path& file)
{
    return file.extension().string() == ".bed";
}

void seek_past_bed_header(std::ifstream& file)
{
    // TODO
}

auto open(const boost::filesystem::path& file)
{
    std::ifstream result {file.string()};
    
    if (is_bed_file(file)) {
        seek_past_bed_header(result);
    }
    
    return result;
}

bool is_valid_bed_record(const std::string& line)
{
    constexpr static char bedDelim {'\t'};
    return std::count(std::cbegin(line), std::cend(line), bedDelim) >= 3;
}

std::string convert_bed_line_to_region_str(const std::string& bed_line)
{
    if (!is_valid_bed_record(bed_line)) {
        throw std::runtime_error {"BadBEDRecord: insufficient columns"};
    }
    
    constexpr static char bedDelim {'\t'};
    
    const auto tokens = utils::split(bed_line, bedDelim);
    
    return std::string {tokens[0] + ':' + tokens[1] + '-' + tokens[2]};
}

std::function<GenomicRegion(const std::string&)>
make_region_line_parser(const boost::filesystem::path& file, const ReferenceGenome& reference)
{
    if (is_bed_file(file)) {
        return [&] (const std::string& line) -> GenomicRegion
        {
            return io::parse_region(convert_bed_line_to_region_str(line), reference);
        };
    } else {
        return [&] (const std::string& line) { return io::parse_region(line, reference); };
    }
}

} // namespace

std::deque<GenomicRegion> extract_regions(const boost::filesystem::path& file,
                                          const ReferenceGenome& reference)
{
    auto stream = open(file);
    
    std::deque<GenomicRegion> result {};
    
    std::transform(std::istream_iterator<Line>(stream), std::istream_iterator<Line>(),
                   std::back_inserter(result), make_region_line_parser(file, reference));
    
    result.shrink_to_fit();
    
    return result;
}

} // namespace io
} // namespace octopus
