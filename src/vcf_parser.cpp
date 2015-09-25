//
//  vcf_parser.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "vcf_parser.hpp"

#include <string>
#include <algorithm> // std::count_if, std::copy, std::transform, std::for_each, std::count
#include <iterator>  // std::cbegin, std::cend, std::advance, std::next
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

VcfHeader VcfParser::fetch_header()
{
    return header_;
}

size_t VcfParser::num_records()
{
    reset_vcf();
    
    return std::count_if(std::istreambuf_iterator<char>(file_), std::istreambuf_iterator<char>(),
                         [] (char c) {
                             return c == '\n';
                         });
}

size_t VcfParser::num_records(const GenomicRegion& region)
{
    reset_vcf();
    
    return std::count_if(std::istream_iterator<Line>(file_), std::istream_iterator<Line>(),
                         [&region] (const auto& line) {
                             return overlaps(line, region);
                         });
}

std::vector<VcfRecord> VcfParser::fetch_records(Unpack level)
{
    reset_vcf();
    
    std::vector<VcfRecord> result {};
    
    bool unpack_all {level == Unpack::All};
    
    std::transform(std::istream_iterator<Line>(file_), std::istream_iterator<Line>(),
                   std::back_inserter(result), [this, unpack_all] (const auto& line) {
                       return (unpack_all) ? parse_record(line, samples_) : parse_record(line);
                   });
    
    return result;
}

std::vector<VcfRecord> VcfParser::fetch_records(const GenomicRegion& region, Unpack level)
{
    reset_vcf();
    
    std::vector<VcfRecord> result {};
    
    bool unpack_all {level == Unpack::All};
    
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
    auto it = std::find(line.cbegin(), line.cend(), '=');
    return it != line.cend() && std::next(it) != line.cend() && *std::next(it) == '<' && line.back() == '>';
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
    if (std::count(field.cbegin(), field.cend(), '=') != 1) {
        throw std::runtime_error {"VCF header field " + field + " is incorrectly formatted"};
    }
    
    auto pos = field.find_first_of('=');
    return std::make_pair(field.substr(0, pos), field.substr(pos + 1));
}

std::unordered_map<std::string, std::string> parse_fields(const std::string& fields)
{
    std::istringstream ss {fields};
    std::unordered_map<std::string, std::string> result {};
    std::transform(std::istream_iterator<Field>(ss), std::istream_iterator<Field>(), std::inserter(result, result.begin()),
                   [] (const auto& field) { return parse_field(field); });
    return result;
}

// ##TAG=<keyA=valueA,...,keyB=valueB>
void parse_structured_header_line(const std::string& line, VcfHeader::Builder& hb)
{
    try {
        auto pos = line.find_first_of('=');
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
    
    std::istream_iterator<Column> eos;
    std::istream_iterator<Column> it {ss};
    
    std::advance(it, 8);
    
    if (it != eos) {
        std::advance(it, 1); // will be on FORMAT
        
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
    
    std::string line {};
    
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

bool overlaps(const std::string& line, const GenomicRegion& region)
{
    std::istringstream ss {line};
    
    std::istream_iterator<Column> it {ss};
    
    if (it->data != region.get_contig_name()) return false; // CHROM
    
    std::advance(it, 1); // POS
    
    auto begin = std::stol(it->data);
    
    std::advance(it, 2); // REF
    
    auto end = begin + static_cast<long>(it->data.length());
    
    return (std::min(static_cast<long>(region.get_end()), end) -
                std::max(static_cast<long>(region.get_begin()), begin)) > 0;
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
        auto pos = f.find_first_of('=');
        rb.add_info(f.substr(0, pos), split(f.substr(pos + 1), ','));
    }
}

void parse_sample(const std::string& column, const VcfRecord::SampleIdType& sample,
                  const std::vector<std::string>& format, VcfRecord::Builder& rb)
{
    auto values = split(column, ':');
    
    if (format.front() == "GT") {
        const std::string& genotype = values.front();
        bool is_phased {(genotype.find_first_of('|') != std::string::npos)};
        auto alleles = split(genotype, (is_phased) ? '|' : '/');
        
        std::vector<unsigned> allele_numbers {};
        allele_numbers.reserve(alleles.size());
        std::transform(alleles.cbegin(), alleles.cend(), std::back_inserter(allele_numbers),
                       [] (const std::string& a) {
                           return static_cast<unsigned>(std::stoul(a));
                       });
        
        rb.add_genotype(sample, allele_numbers, is_phased);
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
    VcfRecord::Builder rb {};
    
    std::istringstream ss {line};
    
    std::istream_iterator<Column> eos;
    std::istream_iterator<Column> it {ss};
    
    rb.set_chromosome(it->data);
    std::advance(it, 1);
    rb.set_position(static_cast<VcfRecord::SizeType>(std::stol(it->data)));
    std::advance(it, 1);
    rb.set_id(it->data);
    std::advance(it, 1);
    rb.set_ref_allele(it->data);
    std::advance(it, 1);
    rb.set_alt_alleles(split(it->data, ','));
    std::advance(it, 1);
    rb.set_quality(static_cast<VcfRecord::QualityType>(std::stoi(it->data)));
    std::advance(it, 1);
    rb.set_filters(split(it->data, ':'));
    std::advance(it, 1);
    parse_info(it->data, rb);
    
    std::advance(it, 1);
    if (!samples.empty() && it != eos) {
        auto format = split(it->data, ':');
        std::advance(it, 1);
        for (const auto& sample : samples) {
            parse_sample(it->data, sample, format, rb); // set after so can move
            std::advance(it, 1);
        }
        rb.set_format(std::move(format));
    }
    
    return rb.build_once();
}
