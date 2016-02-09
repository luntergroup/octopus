//
//  vcf_header.cpp
//  Octopus
//
//  Created by Daniel Cooke on 28/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "vcf_header.hpp"

#include <algorithm>
#include <utility>
#include <iterator>

bool is_valid_line(const std::string& line);
bool is_valid_field(const std::string& str);
bool is_format_line(const std::string& line);
std::unordered_map<std::string, std::string> parse_fields(const std::string& fields);

// public methods

VcfHeader::VcfHeader(std::string file_format)
: file_format_ {std::move(file_format)}
{}

const VcfHeader::ValueType& VcfHeader::get_file_format() const noexcept
{
    return file_format_;
}

unsigned VcfHeader::num_samples() const noexcept
{
    return static_cast<unsigned>(samples_.size());
}

std::vector<std::string> VcfHeader::get_samples() const
{
    return samples_;
}

bool VcfHeader::has_basic_field(const KeyType& key) const noexcept
{
    return basic_fields_.count(key) == 1;
}

bool VcfHeader::has_structured_field(const TagType& tag) const noexcept
{
    return structured_fields_.count(tag) > 0;
}

bool VcfHeader::has_structured_field(const TagType& tag, const KeyType& key) const noexcept
{
    return structured_fields_.count(tag) > 0; // TODO: complete this
}

std::vector<VcfHeader::KeyType> VcfHeader::get_basic_field_keys() const
{
    std::vector<KeyType> result {};
    result.reserve(basic_fields_.size());
    
    std::transform(std::cbegin(basic_fields_), std::cend(basic_fields_),
                   std::back_inserter(result), [] (const auto& p) { return p.first; });
    
    return result;
}

std::vector<VcfHeader::TagType> VcfHeader::get_structured_field_tags() const
{
    std::vector<TagType> result {};
    result.reserve(structured_fields_.size());
    
    std::transform(std::cbegin(structured_fields_), std::cend(structured_fields_),
                   std::back_inserter(result), [] (const auto& p) { return p.first; });
    
    std::sort(std::begin(result), std::end(result));
    
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    
    return result;
}

const VcfHeader::ValueType& VcfHeader::get_basic_field_value(const KeyType& key) const
{
    return basic_fields_.at(key);
}

const VcfHeader::ValueType& VcfHeader::get_structured_field_value(const TagType& tag,
                                                                  const KeyType& id_key,
                                                                  const ValueType& id_value,
                                                                  const KeyType& lookup_key) const
{
    const auto er = structured_fields_.equal_range(tag);
    return std::find_if(er.first, er.second,
                        [&id_key, &id_value, &lookup_key] (const auto& p) {
                            return p.second.at(id_key) == id_value;
                        })->second.at(lookup_key);
}

const std::unordered_map<VcfHeader::KeyType, VcfHeader::ValueType>& VcfHeader::get_basic_fields() const noexcept
{
    return basic_fields_;
}

std::vector<VcfHeader::StructuredField> VcfHeader::get_structured_fields(const TagType& tag) const
{
    std::vector<StructuredField> result {};
    result.reserve(structured_fields_.count(tag));
    
    const auto er = structured_fields_.equal_range(tag);
    
    std::transform(er.first, er.second, std::back_inserter(result),
                   [] (const auto& p) { return p.second; });
    
    return result;
}

const std::unordered_multimap<VcfHeader::TagType, VcfHeader::StructuredField>&
VcfHeader::get_structured_fields() const noexcept
{
    return structured_fields_;
}

// non-member methods

const VcfHeader::ValueType& get_id_field_value(const VcfHeader& header, const VcfHeader::TagType& tag,
                                               const VcfHeader::ValueType& id_value,
                                               const VcfHeader::KeyType& lookup_key)
{
    return header.get_structured_field_value(tag, "ID", id_value, lookup_key);
}

const VcfHeader::ValueType& get_id_field_type(const VcfHeader& header, const VcfHeader::TagType& tag,
                                              const VcfHeader::ValueType& id_value)
{
    return header.get_structured_field_value(tag, "ID", id_value, "Type");
}

VcfType get_typed_value(const VcfHeader& header, const VcfHeader::TagType& tag,
                        const VcfHeader::KeyType& key, const VcfHeader::ValueType& value)
{
    return make_vcf_type(get_id_field_type(header, tag, key), value);
}

VcfType get_typed_info_value(const VcfHeader& header, const VcfHeader::KeyType& key,
                             const VcfHeader::ValueType& value)
{
    return get_typed_value(header, "INFO", key, value);
}

VcfType get_typed_format_value(const VcfHeader& header, const VcfHeader::KeyType& key,
                               const VcfHeader::ValueType& value)
{
    return get_typed_value(header, "FORMAT", key, value);
}

std::vector<VcfType> get_typed_values(const VcfHeader& header, const VcfHeader::KeyType& format_key,
                                      const VcfHeader::KeyType& field_key,
                                      const std::vector<VcfHeader::ValueType>& values)
{
    std::vector<VcfType> result {};
    result.reserve(values.size());
    
    std::transform(values.cbegin(), values.cend(), std::back_inserter(result),
                   [&header, &format_key, &field_key] (const auto& value) {
                       return get_typed_value(header, format_key, field_key, value);
                   });
    
    return result;
}

std::vector<VcfType> get_typed_info_values(const VcfHeader& header, const VcfHeader::KeyType& field_key,
                                           const std::vector<VcfHeader::ValueType>& values)
{
    return get_typed_values(header, "INFO", field_key, values);
}

std::vector<VcfType> get_typed_format_values(const VcfHeader& header, const VcfHeader::KeyType& field_key,
                                             const std::vector<VcfHeader::ValueType>& values)
{
    return get_typed_values(header, "FORMAT", field_key, values);
}

bool contig_line_exists(const VcfHeader& header, const std::string& contig)
{
    return header.has_structured_field("contig", contig);
}

bool operator==(const VcfHeader& lhs, const VcfHeader& rhs)
{
    return lhs.get_file_format() == rhs.get_file_format() && lhs.get_samples() == rhs.get_samples()
            && lhs.get_basic_fields() == rhs.get_basic_fields()
            && lhs.get_structured_fields() == rhs.get_structured_fields();
}

bool operator!=(const VcfHeader& lhs, const VcfHeader& rhs)
{
    return !operator==(lhs, rhs);
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
    os << "##fileformat=" << header.file_format_ << std::endl;
    
    for (const auto& field : header.basic_fields_) {
        os << "##" << field.first << "=" << field.second << std::endl;
    }
    
    for (const auto& format : header.structured_fields_) {
        os << "##" << format.first << "=" << format.second << std::endl;
    }
    
    static const std::vector<std::string> columns {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"};
    
    
    return os;
}

// VcfHeader::Builder

VcfHeader::Builder::Builder(const VcfHeader& header)
:
file_format_ {header.file_format_},
samples_ {header.samples_},
basic_fields_ {header.basic_fields_},
structured_fields_ {header.structured_fields_}
{}

VcfHeader::Builder& VcfHeader::Builder::set_file_format(std::string file_format)
{
    file_format_ = std::move(file_format);
    return *this;
}

VcfHeader::Builder& VcfHeader::Builder::add_sample(std::string sample)
{
    samples_.push_back(std::move(sample));
    return *this;
}

VcfHeader::Builder& VcfHeader::Builder::set_samples(std::vector<std::string> samples)
{
    samples_ = std::move(samples);
    return *this;
}

VcfHeader::Builder& VcfHeader::Builder::add_basic_field(std::string key, std::string value)
{
    if (key != "fileformat") {
        basic_fields_.emplace(std::move(key), std::move(value));
    }
    return *this;
}

VcfHeader::Builder& VcfHeader::Builder::add_structured_field(std::string tag,
                                                             std::unordered_map<std::string, std::string> values)
{
    structured_fields_.emplace(std::move(tag), std::move(values));
    return *this;
}

std::string add_quotes(const std::string& str)
{
    std::string result {};
    if (str.front() != '"') result = '"' + str;
    if (str.back() != '"') result.push_back('"');
    return result;
}

VcfHeader::Builder& VcfHeader::Builder::add_info(std::string id, std::string number,
                                                 std::string type, std::string description,
                                                 std::unordered_map<std::string, std::string> other_values)
{
    other_values.reserve(other_values.size() + 4);
    
    other_values.emplace("ID", std::move(id));
    other_values.emplace("Number", std::move(number));
    other_values.emplace("Type", std::move(type));
    other_values.emplace("Description", add_quotes(description));
    
    structured_fields_.emplace("INFO", std::move(other_values));
    
    return *this;
}

VcfHeader::Builder& VcfHeader::Builder::add_filter(std::string id, std::string description,
                                                   std::unordered_map<std::string, std::string> other_values)
{
    other_values.reserve(other_values.size() + 2);
    
    other_values.emplace("ID", std::move(id));
    other_values.emplace("Description", add_quotes(description));
    
    structured_fields_.emplace("FILTER", std::move(other_values));
    
    return *this;
}

VcfHeader::Builder& VcfHeader::Builder::add_format(std::string id, std::string number,
                                                   std::string type, std::string description,
                                                   std::unordered_map<std::string, std::string> other_values)
{
    other_values.reserve(other_values.size() + 4);
    
    other_values.emplace("ID", std::move(id));
    other_values.emplace("Number", std::move(number));
    other_values.emplace("Type", std::move(type));
    other_values.emplace("Description", add_quotes(description));
    
    structured_fields_.emplace("FORMAT", std::move(other_values));
    
    return *this;
}

VcfHeader::Builder& VcfHeader::Builder::add_contig(std::string id,
                                                   std::unordered_map<std::string, std::string> other_values)
{
    other_values.emplace("ID", std::move(id));
    
    structured_fields_.emplace("contig", std::move(other_values));
    
    return *this;
}

VcfHeader VcfHeader::Builder::build() const
{
    return VcfHeader {file_format_, samples_, basic_fields_, structured_fields_};
}

VcfHeader VcfHeader::Builder::build_once() noexcept
{
    return VcfHeader {std::move(file_format_), std::move(samples_), std::move(basic_fields_), std::move(structured_fields_)};
}

VcfHeader::Builder get_default_header_builder()
{
    VcfHeader::Builder result {};
    
    result.add_info("AA", "1", "String", "Ancestral allele");
    result.add_info("AC", "1", "Integer", "Allele count in genotypes, for each ALT allele, in the same order as listed");
    result.add_info("AF", "A", "Float", "Allele Frequency, for each ALT allele, in the same order as listed");
    result.add_info("AN", "1", "Integer", "Total number of alleles in called genotypes");
    result.add_info("BQ", "1", "Integer", "RMS base quality at this position");
    result.add_info("CIGAR", "A", "String", "Cigar string describing how to align an alternate allele to the reference allele");
    result.add_info("DB", "0", "Flag", "dbSNP membership");
    result.add_info("DP", "1", "Integer", "Combined depth across samples");
    result.add_info("END", "1", "Integer", "End position of the variant described in this record");
    result.add_info("H2", "0", "Flag", "Membership in hapmap2");
    result.add_info("H3", "0", "Flag", "Membership in hapmap3");
    result.add_info("MQ", "1", "Integer", "RMS mapping quality");
    result.add_info("MQ0", "1", "Integer", "Number of MAPQ == 0 reads covering this record");
    result.add_info("NS", "1", "Integer", "Number of samples with data");
    result.add_info("SB", "1", "Integer", "Strand bias at this position");
    result.add_info("SOMATIC", "0", "Flag", "Indicates that the record is a somatic mutation, for cancer genomics");
    result.add_info("VALIDATED", "0", "Flag", "Validated by follow-up experiment");
    result.add_info("1000G", "0", "Flag", "Membership in 1000 Genomes");
    
    result.add_format("GT", "1", "String", "Genotype");
    result.add_format("DP", "1", "Integer", "Read depth at this position for this sample");
    result.add_format("FT", "1", "String", "Sample genotype filter indicating if this genotype was “called”");
    result.add_format("GL", "G", "Float", "log10-scaled genotype likelihoods");
    result.add_format("GLE", "1", "Integer", "Genotype likelihoods of heterogeneous ploidy");
    result.add_format("PL", "G", "Integer", "Phred-scaled genotype likelihoods");
    result.add_format("GP", "G", "Float", "Phred-scaled genotype posterior probabilities");
    result.add_format("GQ", "1", "Integer", "Conditional genotype quality (phred-scaled)");
    result.add_format("HQ", "1", "Integer", "Haplotype qualities");
    result.add_format("PS", "1", "String", "Phase set");
    result.add_format("PQ", "1", "Integer", "Phasing quality");
    result.add_format("EC", "1", "Integer", "Expected alternate allele counts");
    result.add_format("MQ", "1", "Integer", "RMS mapping quality");
    result.add_format("BQ", "1", "Integer", "RMS base quality at this position");
    
    result.add_filter("PASS", "All filters passed");
    result.add_filter("REFCALL", "All samples are called homozygous reference at the site/block");
    
    return result;
}
