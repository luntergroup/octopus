//
//  vcf_parser.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "vcf_parser.hpp"

#include <algorithm>
#include <iterator>
#include <stdexcept>

#include <boost/optional.hpp>
#include <boost/lexical_cast.hpp>

#include "genomic_region.hpp"
#include "vcf_record.hpp"

#include <iostream> // TEST

VcfHeader parse_header(std::ifstream& vcf_file);
bool overlaps(const std::string& line, const GenomicRegion& region);
VcfRecord parse_record(const std::string& line, const std::vector<VcfRecord::SampleIdType>& samples = {});

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

using Line   = Token<'\n'>;
using Column = Token<'\t'>;

// public methods

VcfParser::VcfParser(const fs::path& file_path)
:
file_path_ {file_path},
file_ {file_path_.string()},
header_ {parse_header(file_)},
samples_ {header_.samples()},
first_record_pos_ {file_.tellg()}
{}

bool VcfParser::is_header_written() const noexcept
{
    return true; // always the case as can only read
}

VcfHeader VcfParser::fetch_header() const
{
    return header_;
}

std::size_t VcfParser::count_records() const
{
    const auto result = std::count_if(std::istreambuf_iterator<char>(file_), std::istreambuf_iterator<char>(),
                                      [] (char c) { return c == '\n'; });
    
    reset_vcf();
    
    return result;
}

std::size_t VcfParser::count_records(const std::string& contig) const
{
    const auto result = std::count_if(std::istream_iterator<Line>(file_), std::istream_iterator<Line>(),
                                      [&contig] (const auto& line) { return is_same_contig(line, contig); });
    
    reset_vcf();
    
    return result;
}

std::size_t VcfParser::count_records(const GenomicRegion& region) const
{
    const auto result = std::count_if(std::istream_iterator<Line>(file_), std::istream_iterator<Line>(),
                                      [&region] (const auto& line) { return overlaps(line, region); });
    
    reset_vcf();
    
    return result;
}

std::vector<VcfRecord> VcfParser::fetch_records(const UnpackPolicy level) const
{
    std::vector<VcfRecord> result {};
    
    result.reserve(count_records());
    
    bool unpack_all {level == UnpackPolicy::All};
    
    std::transform(std::istream_iterator<Line>(file_), std::istream_iterator<Line>(),
                   std::back_inserter(result), [this, unpack_all] (const auto& line) {
                       return (unpack_all) ? parse_record(line, samples_) : parse_record(line);
                   });
    
    reset_vcf();
    
    return result;
}

std::vector<VcfRecord> VcfParser::fetch_records(const std::string& contig, const UnpackPolicy level) const
{
    std::vector<VcfRecord> result {};
    
    result.reserve(count_records(contig));
    
    bool unpack_all {level == UnpackPolicy::All};
    
    std::for_each(std::istream_iterator<Line>(file_), std::istream_iterator<Line>(),
                  [this, &result, &contig, unpack_all] (const auto& line) {
                      if (is_same_contig(line, contig)) {
                          result.push_back((unpack_all) ? parse_record(line, samples_) : parse_record(line));
                      }
                  });
    
    reset_vcf();
    
    return result;
}

std::vector<VcfRecord> VcfParser::fetch_records(const GenomicRegion& region, const UnpackPolicy level) const
{
    std::vector<VcfRecord> result {};
    
    result.reserve(count_records(region));
    
    bool unpack_all {level == UnpackPolicy::All};
    
    std::for_each(std::istream_iterator<Line>(file_), std::istream_iterator<Line>(),
                  [this, &result, &region, unpack_all] (const std::string& line) {
                      if (overlaps(line, region)) {
                          result.push_back((unpack_all) ? parse_record(line, samples_) : parse_record(line));
                      }
                  });
    
    reset_vcf();
    
    return result;
}

// private methods

void VcfParser::reset_vcf() const
{
    file_.clear();
    file_.seekg(first_record_pos_);
}

// non-member methods

using Field = Token<','>;

// A field looks like key=value, fields are delimited with ',' e.g. keyA=valueA,...,keyB=valueB.
// But value may be quoted (i.e. "value"), and anything goes inside the quotes. So we must make
// sure we actually find the next ',' for the next Field, rather than a ',' inside the quotes.
std::istream& operator>>(std::istream& str, Field& field)
{
    std::getline(str, field.data, ',');
    
    auto pos = field.data.find_first_of('=');
    
    if (pos != field.data.length() - 1 && field.data[pos + 1] == '"') {
        std::string s;
        while (field.data.back() != '"') {
            std::getline(str, s);
            field.data += s;
        }
    }
    
    return str;
}

bool is_header_meta_line(const std::string& line)
{
    return line.length() > 3 && line[0] == '#' && line[1] == '#';
}

bool is_structured_header_line(const std::string& line)
{
    auto it = std::find(std::cbegin(line), std::cend(line), '=');
    return it != std::cend(line) && std::next(it) != std::cend(line)
            && *std::next(it) == '<' && line.back() == '>';
}

// ##key=value
void parse_basic_header_line(const std::string& line, VcfHeader::Builder& hb)
{
    if (std::count(line.cbegin(), line.cend(), '=') != 1) {
        throw std::runtime_error {"VCF header line " + line + " is incorrectly formatted"};
    }
    
    auto pos = line.find_first_of('=');
    hb.add_basic_field(line.substr(2, pos - 2), line.substr(pos + 1));
}

std::pair<std::string, std::string> parse_field(const std::string& field)
{
    if (std::count(field.cbegin(), field.cend(), '=') == 0) {
        throw std::runtime_error {"VCF header field " + field + " is incorrectly formatted"};
    }
    auto pos = field.find_first_of('=');
    return std::make_pair(field.substr(0, pos), field.substr(pos + 1));
}

std::unordered_map<std::string, std::string> parse_fields(const std::string& fields)
{
    std::istringstream ss {fields};
    std::unordered_map<std::string, std::string> result {};
    std::transform(std::istream_iterator<Field>(ss), std::istream_iterator<Field>(),
                   std::inserter(result, result.begin()),
                   [] (const auto& field) { return parse_field(field); });
    return result;
}

// ##TAG=<keyA=valueA,...,keyB=valueB>
void parse_structured_header_line(const std::string& line, VcfHeader::Builder& hb)
{
    try {
        const auto pos = line.find_first_of('=');
        auto tag = line.substr(2, pos - 2);
        hb.add_structured_field(std::move(tag), parse_fields(line.substr(pos + 2, line.length() - pos - 3)));
    } catch (...) {
        throw std::runtime_error {"VCF header line " + line + " is incorrectly formatted"};
    }
}

void parse_header_meta_line(const std::string& line, VcfHeader::Builder& hb)
{
    if (is_structured_header_line(line)) {
        parse_structured_header_line(line, hb);
    } else {
        parse_basic_header_line(line, hb);
    }
}

// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\t...\tSAMPLEN
void parse_header_sample_names(const std::string& line, VcfHeader::Builder& hb)
{
    std::istringstream ss {line};
    
    std::istream_iterator<Column> it {ss}, eos {};
    
    std::advance(it, 8);
    
    if (it != eos) {
        std::advance(it, 1); // now on FORMAT
        
        std::vector<std::string> samples {};
        std::copy(it, eos, std::back_inserter(samples));
        samples.shrink_to_fit();
        
        hb.set_samples(std::move(samples));
    }
}

VcfHeader parse_header(std::ifstream& vcf_file)
{
    vcf_file.seekg(0, std::ios::beg); // reset
    
    VcfHeader::Builder hb {};
    
    std::string line;
    std::getline(vcf_file, line);
    
    if (!is_header_meta_line(line)) {
        throw std::runtime_error {"the first line of a VCF file must be ##fileformat"};
    }
    
    hb.set_file_format(line.substr(line.find_first_of('=') + 1));
    
    while (std::getline(vcf_file, line) && is_header_meta_line(line)) {
        parse_header_meta_line(line, hb);
    }
    
    parse_header_sample_names(line, hb); // last line is column names, including sample names
    
    return hb.build_once();
}

bool is_same_contig(const std::string& line, const std::string& contig)
{
    std::istringstream ss {line};
    
    std::istream_iterator<Column> it {ss};
    
    return it->data == contig;
}

bool overlaps(const std::string& line, const GenomicRegion& region)
{
    std::istringstream ss {line};
    
    std::istream_iterator<Column> it {ss};
    
    if (it->data != region.contig_name()) return false; // CHROM
    
    std::advance(it, 1); // POS
    
    const auto begin = std::stol(it->data);
    
    std::advance(it, 2); // REF
    
    const auto end = begin + static_cast<long>(it->data.length());
    
    return (std::min(static_cast<long>(region.end()), end) -
                std::max(static_cast<long>(region.begin()), begin)) > 0;
}

std::vector<std::string> split(const std::string& str, char delim = ',')
{
    std::stringstream ss {str};
    std::string item;
    
    std::vector<std::string> result {};
    result.reserve(std::count(std::cbegin(str), std::cend(str), delim) + 1);
    
    while (std::getline(ss, item, delim)) {
        result.emplace_back(item);
    }
    
    return result;
}

void parse_info_field(const std::string& field, VcfRecord::Builder& rb)
{
    const auto pos = field.find_first_of('=');
    
    if (pos == std::string::npos) {
        rb.add_info(field);
    } else {
        rb.add_info(field.substr(0, pos), split(field.substr(pos + 1), ','));
    }
}

using InfoField = Token<';'>;

void parse_info(const std::string& column, VcfRecord::Builder& rb)
{
    std::istringstream ss {column};
    
    std::for_each(std::istream_iterator<InfoField>(ss), std::istream_iterator<InfoField>(),
                  [&rb] (const auto& field) {
                      parse_info_field(field, rb);
                  });
}

bool is_phased(const std::string& genotype)
{
    const auto pos = genotype.find_first_of("|/");
    
    if (pos == std::string::npos) {
        return true; // must be haploid
    }
    
    return genotype[pos] == '|';
}

void parse_genotype(const VcfRecord::SampleIdType& sample, const std::string& genotype,
                    VcfRecord::Builder& rb)
{
    const bool phased {is_phased(genotype)};
    
    auto allele_numbers = split(genotype, (phased) ? '|' : '/');
    
    using Phasing = VcfRecord::Builder::Phasing;
    
    std::vector<boost::optional<unsigned>> alleles {};
    
    const auto ploidy = allele_numbers.size();
    
    alleles.reserve(ploidy);
    
    std::transform(std::cbegin(allele_numbers), std::cend(allele_numbers), std::back_inserter(alleles),
                   [ploidy] (const auto& a) -> boost::optional<unsigned> {
                       try {
                           const auto i = boost::lexical_cast<unsigned>(a);
                           
                           if (i < ploidy) {
                               return i;
                           } else {
                               return boost::none;
                           }
                       } catch (const boost::bad_lexical_cast&) {
                           return boost::none;
                       }
                   });
    
    rb.add_genotype(sample, std::move(alleles), (phased) ? Phasing::Phased : Phasing::Unphased);
}

using SampleField = Token<':'>;

void parse_sample(const std::string& column, const VcfRecord::SampleIdType& sample,
                  const std::vector<std::string>& format, VcfRecord::Builder& rb)
{
    auto first_key = std::cbegin(format);
    
    std::istringstream ss {column};
    
    std::istream_iterator<SampleField> first_value {ss};
    
    if (format.front() == "GT") { // GT must always come first, if present
        parse_genotype(sample, *first_value, rb);
        
        ++first_key;
        ++first_value;
    }
    
    std::for_each(first_value, std::istream_iterator<SampleField> {},
                  [&rb, &sample, &first_key] (const std::string& value) {
                      rb.add_genotype_field(sample, *first_key, split(value, ','));
                      ++first_key;
                  });
}

VcfRecord parse_record(const std::string& line, const std::vector<VcfRecord::SampleIdType>& samples)
{
    std::istringstream ss {line};
    
    std::istream_iterator<Column> it {ss}, eos {};
    
    VcfRecord::Builder rb {};
    
    rb.set_chrom(it->data);
    ++it;
    rb.set_pos(static_cast<VcfRecord::SizeType>(std::stol(it->data)));
    ++it;
    rb.set_id(it->data);
    ++it;
    rb.set_ref(it->data);
    ++it;
    rb.set_alt(split(it->data, ','));
    ++it;
    
    if (it->data == ".") {
        rb.set_qual(0);
    } else {
        try {
            rb.set_qual(static_cast<VcfRecord::QualityType>(std::stod(it->data)));
        } catch (const std::invalid_argument& e) {
            rb.set_qual(0); // or should throw?
        }
    }
    
    ++it;
    if (it->data != ".") {
        rb.set_filter(split(it->data, ':'));
    }
    ++it;
    parse_info(it->data, rb);
    ++it;
    
    if (!samples.empty() && it != eos) {
        auto format = split(it->data, ':');
        ++it;
        for (const auto& sample : samples) {
            parse_sample(it->data, sample, format, rb); // set after so can move
            ++it;
        }
        rb.set_format(std::move(format));
    }
    
    return rb.build_once();
}
