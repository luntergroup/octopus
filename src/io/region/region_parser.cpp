// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "region_parser.hpp"

#include <sstream>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <functional>

#include <boost/optional.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/type_traits/is_unsigned.hpp>

#include "exceptions/user_error.hpp"
#include "utils/string_utils.hpp"

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
    : region_ {std::move(region)}
    , why_ {std::move(why)}
    {}
    
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
    : region_ {std::move(region)}
    , reference_ {std::move(reference)}
    {}
};

template <bool is_unsigned>
struct unsigned_checker
{
    template <typename String_type>
    static inline void do_check(const String_type& str) {}
};

template <>
struct unsigned_checker<true>
{
    template <typename String_type>
    static inline void do_check(const String_type& str)
    {
        if (str.front() == '-') boost::throw_exception(boost::bad_lexical_cast());
    }
};

GenomicRegion parse_region(std::string region, const ReferenceGenome& reference)
{
    if (reference.has_contig(region)) {
        // contig
        return GenomicRegion {std::move(region), 0, reference.contig_size(region)};
    }
    const auto last_colon_ritr = std::find(std::rbegin(region), std::rend(region), ':');
    if (last_colon_ritr == std::rend(region)) {
        throw MissingReferenceContig {region, reference.name()};
    }
    const auto begin_begin_itr = last_colon_ritr.base();
    region.erase(std::remove(begin_begin_itr, std::end(region), ','), std::end(region));
    const auto begin_end_itr = std::find(begin_begin_itr, std::end(region), '-');
    if (begin_end_itr == begin_begin_itr) {
        throw MalformedRegion {region};
    }
    const std::string begin_str {begin_begin_itr, begin_end_itr};
    try {
        using Position = GenomicRegion::Position;
        std::string contig {std::begin(region), std::prev(begin_begin_itr)};
        if (!reference.has_contig(contig)) {
            throw MissingReferenceContig {region, reference.name()};
        }
        const auto contig_size = static_cast<Position>(reference.contig_size(contig));
        if (contig_size == 0) {
            return GenomicRegion {std::move(contig), 0, 0};
        }
        unsigned_checker<boost::is_unsigned<Position>::value>::do_check(begin_str);
        const auto begin = std::min(boost::lexical_cast<Position>(begin_str), contig_size - 1);
        if (begin_end_itr == std::end(region)) {
            // contig:position
            return GenomicRegion {std::move(contig), begin, begin + 1};
        } else if (std::next(begin_end_itr) == std::end(region)) {
            // contig:begin-
            region.erase(std::prev(begin_begin_itr), std::end(region));
            return GenomicRegion {std::move(region), begin, contig_size};
        } else {
            // contig:begin-end
            const std::string end_str {std::next(begin_end_itr), std::end(region)};
            unsigned_checker<boost::is_unsigned<Position>::value>::do_check(end_str);
            const auto end = std::min(boost::lexical_cast<Position>(end_str), contig_size);
            if (begin > end) {
                throw MalformedRegion {region, "has begin greater than end"};
            }
            return GenomicRegion {std::move(contig), begin, end};
        }
    } catch (const boost::bad_lexical_cast&) {
        throw MalformedRegion {region};
    }
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
    return std::count(std::cbegin(line), std::cend(line), bedDelim) >= 2;
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
make_throwing_region_line_parser(const boost::filesystem::path& file, const ReferenceGenome& reference)
{
    if (is_bed_file(file)) {
        return [&] (const std::string& line) -> GenomicRegion
        {
            return parse_region(convert_bed_line_to_region_str(line), reference);
        };
    } else {
        return [&] (const std::string& line) -> GenomicRegion
        {
            return parse_region(line, reference);
        };
    }
}

std::function<boost::optional<GenomicRegion>(const std::string&)>
make_optional_region_line_parser(const boost::filesystem::path& file, const ReferenceGenome& reference)
{
    if (is_bed_file(file)) {
        return [&] (const std::string& line) -> boost::optional<GenomicRegion>
        {
            try {
                return parse_region(convert_bed_line_to_region_str(line), reference);
            } catch (const MissingReferenceContig&) {
                return boost::none;
            }
        };
    } else {
        return [&] (const std::string& line) -> boost::optional<GenomicRegion>
        {
            try {
                return parse_region(line, reference);
            } catch (const MissingReferenceContig&) {
                return boost::none;
            }
        };
    }
}

} // namespace

std::deque<GenomicRegion> extract_regions(const boost::filesystem::path& file,
                                          const ReferenceGenome& reference,
                                          const NonreferenceContigPolicy policy)
{
    auto stream = open(file);
    std::deque<GenomicRegion> result {};
    if (policy == NonreferenceContigPolicy::exception) {
        std::transform(std::istream_iterator<Line>(stream), std::istream_iterator<Line>(),
                       std::back_inserter(result), make_throwing_region_line_parser(file, reference));
    } else {
        std::deque<boost::optional<GenomicRegion>> filtered_results {};
        std::transform(std::istream_iterator<Line>(stream), std::istream_iterator<Line>(),
                       std::back_inserter(filtered_results), make_optional_region_line_parser(file, reference));
        for (const auto& region : filtered_results) {
            if (region) {
                result.push_back(std::move(*region));
            }
        }
    }
    result.shrink_to_fit();
    return result;
}

} // namespace io
} // namespace octopus
