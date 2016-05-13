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
samples_ {header_.get_samples()},
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

std::size_t VcfParser::count_records()
{
    reset_vcf();
    
    return std::count_if(std::istreambuf_iterator<char>(file_), std::istreambuf_iterator<char>(),
                         [] (char c) { return c == '\n'; });
}

std::size_t VcfParser::count_records(const std::string& contig)
{
    reset_vcf();
    
    return std::count_if(std::istream_iterator<Line>(file_), std::istream_iterator<Line>(),
                         [&contig] (const auto& line) { return is_same_contig(line, contig); });
}

std::size_t VcfParser::count_records(const GenomicRegion& region)
{
    reset_vcf();
    
    return std::count_if(std::istream_iterator<Line>(file_), std::istream_iterator<Line>(),
                         [&region] (const auto& line) { return overlaps(line, region); });
}

std::vector<VcfRecord> VcfParser::fetch_records(const UnpackPolicy level)
{
    reset_vcf();
    
    std::vector<VcfRecord> result {};
    
    bool unpack_all {level == UnpackPolicy::All};
    
    std::transform(std::istream_iterator<Line>(file_), std::istream_iterator<Line>(),
                   std::back_inserter(result), [this, unpack_all] (const auto& line) {
                       return (unpack_all) ? parse_record(line, samples_) : parse_record(line);
                   });
    
    return result;
}

std::vector<VcfRecord> VcfParser::fetch_records(const std::string& contig, const UnpackPolicy level)
{
    reset_vcf();
    
    std::vector<VcfRecord> result {};
    
    bool unpack_all {level == UnpackPolicy::All};
    
    std::for_each(std::istream_iterator<Line>(file_), std::istream_iterator<Line>(),
                  [this, &result, &contig, unpack_all] (const auto& line) {
                      if (is_same_contig(line, contig)) {
                          result.emplace_back((unpack_all) ? parse_record(line, samples_) : parse_record(line));
                      }
                  });
    
    return result;
}

std::vector<VcfRecord> VcfParser::fetch_records(const GenomicRegion& region, const UnpackPolicy level)
{
    reset_vcf();
    
    std::vector<VcfRecord> result {};
    
    bool unpack_all {level == UnpackPolicy::All};
    
    std::for_each(std::istream_iterator<Line>(file_), std::istream_iterator<Line>(),
                  [this, &result, &region, unpack_all] (const std::string& line) {
                      if (overlaps(line, region)) {
                          result.emplace_back((unpack_all) ? parse_record(line, samples_) : parse_record(line));
                      }
                  });
    
    return result;
}

// private methods

void VcfParser::reset_vcf()
{
    file_.clear();
    file_.seekg(first_record_pos_);
}

// non-member methods

using Field  = Token<','>;

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
    std::vector<std::string> result {};
    result.reserve(std::count(str.cbegin(), str.cend(), delim));
    std::stringstream ss {str};
    std::string item;
    
    while (std::getline(ss, item, delim)) {
        result.emplace_back(item);
    }
    
    return result;
}

void parse_info(const std::string& column, VcfRecord::Builder& rb)
{
    for (auto& f : split(column, ';')) {
        const auto pos = f.find_first_of('=');
        rb.add_info(f.substr(0, pos), split(f.substr(pos + 1), ','));
    }
}

void parse_sample(const std::string& column, const VcfRecord::SampleIdType& sample,
                  const std::vector<std::string>& format, VcfRecord::Builder& rb)
{
    auto values = split(column, ':');
    
    if (format.front() == "GT") {
        const std::string& genotype = values.front();
        const bool is_phased {(genotype.find_first_of('|') != std::string::npos)};
        
        auto alleles = split(genotype, (is_phased) ? '|' : '/');
        
        std::vector<unsigned> allele_numbers {};
        allele_numbers.reserve(alleles.size());
        std::transform(std::cbegin(alleles), std::cend(alleles), std::back_inserter(allele_numbers),
                       [] (const std::string& a) { return static_cast<unsigned>(std::stoul(a)); });
        
        rb.add_genotype(sample, allele_numbers,
                        (is_phased) ? VcfRecord::Builder::Phasing::Phased : VcfRecord::Builder::Phasing::Unphased);
    }
    
    auto key_it = std::cbegin(format);
    std::for_each((format.front() == "GT") ? std::cbegin(values) : std::next(std::cbegin(values)), std::cend(values),
                  [&rb, &sample, &key_it] (const std::string& value) {
                      rb.add_genotype_field(sample, *key_it, split(value, ','));
                      ++key_it;
                  });
}

VcfRecord parse_record(const std::string& line, const std::vector<VcfRecord::SampleIdType>& samples)
{
    std::istringstream ss {line};
    
    std::istream_iterator<Column> it {ss}, eos {};
    
    VcfRecord::Builder rb {};
    
    rb.set_chromosome(it->data);
    ++it;
    rb.set_position(static_cast<VcfRecord::SizeType>(std::stol(it->data)));
    ++it;
    rb.set_id(it->data);
    ++it;
    rb.set_ref_allele(it->data);
    ++it;
    rb.set_alt_alleles(split(it->data, ','));
    ++it;
    
    if (it->data == ".") {
        rb.set_quality(0);
    } else {
        try {
            rb.set_quality(static_cast<VcfRecord::QualityType>(std::stoi(it->data)));
        } catch (const std::invalid_argument& e) {
            rb.set_quality(0); // or should throw?
        }
    }
    
    ++it;
    rb.set_filters(split(it->data, ':'));
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
