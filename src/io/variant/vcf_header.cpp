// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "vcf_header.hpp"

#include <algorithm>
#include <utility>
#include <iterator>

#include "vcf_spec.hpp"

namespace octopus {

bool is_valid_line(const std::string& line);
bool is_valid_field(const std::string& str);
bool is_format_line(const std::string& line);
std::unordered_map<std::string, std::string> parse_fields(const std::string& fields);

// public methods

VcfHeader::VcfHeader() : VcfHeader {vcfspec::version} {}

VcfHeader::VcfHeader(std::string file_format)
: file_format_ {std::move(file_format)}
{}

const VcfHeader::Value& VcfHeader::file_format() const noexcept
{
    return file_format_;
}

unsigned VcfHeader::num_samples() const noexcept
{
    return static_cast<unsigned>(samples_.size());
}

const std::vector<std::string>& VcfHeader::samples() const noexcept
{
    return samples_;
}

bool VcfHeader::has(const BasicKey& k) const noexcept
{
    return basic_fields_.count(k) == 1;
}

bool VcfHeader::has(const Tag& t) const noexcept
{
    return structured_fields_.count(t) > 0;
}

bool VcfHeader::has(const Tag& t, const StructuredKey& k) const noexcept
{
    const auto er = structured_fields_.equal_range(t);
    return std::find_if(er.first, er.second,
                        [&] (const auto& p) {
                            return p.second.at(vcfspec::header::meta::struc::id) == k;
                        }) != er.second;
}

const VcfHeader::Value& VcfHeader::at(const BasicKey& k) const
{
    return basic_fields_.at(k);
}

const VcfHeader::Value& VcfHeader::find(const Tag& search_tag, const StructuredKey& id_key,
                                        const Value& id_value, const StructuredKey& search_key) const
{
    const auto er = structured_fields_.equal_range(search_tag);
    return std::find_if(er.first, er.second,
                        [&] (const auto& p) {
                            return p.second.at(id_key) == id_value;
                        })->second.at(search_key);
}

std::vector<VcfHeader::BasicKey> VcfHeader::basic_keys() const
{
    std::vector<BasicKey> result {};
    result.reserve(basic_fields_.size());
    std::transform(std::cbegin(basic_fields_), std::cend(basic_fields_),
                   std::back_inserter(result), [] (const auto& p) { return p.first; });
    return result;
}

std::vector<VcfHeader::Tag> VcfHeader::tags() const
{
    std::vector<Tag> result {};
    result.reserve(structured_fields_.size());
    std::transform(std::cbegin(structured_fields_), std::cend(structured_fields_),
                   std::back_inserter(result), [] (const auto& p) { return p.first; });
    std::sort(std::begin(result), std::end(result));
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

std::vector<VcfHeader::StructuredKey> VcfHeader::keys(const Tag& t) const
{
    std::vector<StructuredKey> result {};
    // TODO
    return result;
}

const VcfHeader::BasicFieldMap& VcfHeader::basic_fields() const noexcept
{
    return basic_fields_;
}

std::vector<VcfHeader::StructuredField> VcfHeader::structured_fields(const Tag& tag) const
{
    std::vector<StructuredField> result {};
    result.reserve(structured_fields_.count(tag));
    const auto er = structured_fields_.equal_range(tag);
    std::transform(er.first, er.second, std::back_inserter(result),
                   [] (const auto& p) { return p.second; });
    return result;
}

const VcfHeader::StructuredFieldMap& VcfHeader::structured_fields() const noexcept
{
    return structured_fields_;
}

// non-member methods

const VcfHeader::Value&
get_id_field_value(const VcfHeader& header, 
                   const VcfHeader::Tag& tag,
                   const VcfHeader::StructuredKey& id_value,
                   const VcfHeader::StructuredKey& lookup_key)
{
    return header.find(tag, vcfspec::header::meta::struc::id, id_value, lookup_key);
}

const VcfHeader::Value& 
get_id_field_type(const VcfHeader& header, 
                  const VcfHeader::Tag& tag,
                  const VcfHeader::Value& id_value)
{
    using namespace vcfspec::header::meta;
    return header.find(tag, struc::id, id_value, struc::type);
}

VcfType get_typed_value(const VcfHeader& header, const VcfHeader::Tag& tag,
                        const VcfHeader::StructuredKey& key, const VcfHeader::Value& value)
{
    return make_vcf_type(get_id_field_type(header, tag, key), value);
}

VcfType get_typed_info_value(const VcfHeader& header, const VcfHeader::StructuredKey& key,
                             const VcfHeader::Value& value)
{
    return get_typed_value(header, vcfspec::header::meta::tag::info, key, value);
}

VcfType get_typed_format_value(const VcfHeader& header, const VcfHeader::StructuredKey& key,
                               const VcfHeader::Value& value)
{
    return get_typed_value(header, vcfspec::header::meta::tag::format, key, value);
}

std::vector<VcfType> get_typed_values(const VcfHeader& header, const VcfHeader::StructuredKey& format_key,
                                      const VcfHeader::StructuredKey& field_key,
                                      const std::vector<VcfHeader::Value>& values)
{
    std::vector<VcfType> result {};
    result.reserve(values.size());
    std::transform(values.cbegin(), values.cend(), std::back_inserter(result),
                   [&header, &format_key, &field_key] (const auto& value) {
                       return get_typed_value(header, format_key, field_key, value);
                   });
    return result;
}

std::vector<VcfType> get_typed_info_values(const VcfHeader& header, const VcfHeader::StructuredKey& field_key,
                                           const std::vector<VcfHeader::Value>& values)
{
    return get_typed_values(header, vcfspec::header::meta::tag::info, field_key, values);
}

std::vector<VcfType> get_typed_format_values(const VcfHeader& header, const VcfHeader::StructuredKey& field_key,
                                             const std::vector<VcfHeader::Value>& values)
{
    return get_typed_values(header, vcfspec::header::meta::tag::format, field_key, values);
}

bool contig_line_exists(const VcfHeader& header, const std::string& contig)
{
    return header.has(vcfspec::header::meta::tag::contig, contig);
}

bool operator==(const VcfHeader& lhs, const VcfHeader& rhs)
{
    return lhs.file_format() == rhs.file_format()
            && lhs.samples() == rhs.samples()
            && lhs.basic_fields() == rhs.basic_fields()
            && lhs.structured_fields() == rhs.structured_fields();
}

bool operator!=(const VcfHeader& lhs, const VcfHeader& rhs)
{
    return !operator==(lhs, rhs);
}

std::ostream& operator<<(std::ostream& os, const VcfHeader::StructuredField& fields)
{
    os << "<";
    for (const auto& p : fields) {
        os << p.first.value << "=" << p.second << ",";
    }
    os << ">";
    return os;
}

std::ostream& operator<<(std::ostream& os, const VcfHeader& header)
{
    using vcfspec::header::lineOpener;
    os << lineOpener << vcfspec::header::meta::vcfVersion;
    os << header.file_format_ << std::endl;
    for (const auto& field : header.basic_fields_) {
        os << lineOpener << field.first.value << "=" << field.second << std::endl;
    }
    for (const auto& format : header.structured_fields_) {
        os << lineOpener << format.first.value << "=" << format.second << std::endl;
    }
    // TODO
    return os;
}

// VcfHeader::Builder

VcfHeader::Builder::Builder(const VcfHeader& header)
: file_format_ {header.file_format_}
, samples_ {header.samples_}
, basic_fields_ {header.basic_fields_}
, structured_fields_ {header.structured_fields_}
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
    if (key != vcfspec::header::meta::vcfVersion) {
        basic_fields_.emplace(std::move(key), std::move(value));
    }
    return *this;
}

VcfHeader::Builder& VcfHeader::Builder::add_structured_field(std::string tag,
                                                             std::unordered_map<std::string, std::string> values)
{
    VcfHeader::StructuredField tmp {};
    tmp.reserve(values.size());
    for (auto&& p : values) tmp.emplace(p.first, std::move(p.second));
    structured_fields_.emplace(std::move(tag), std::move(tmp));
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
    using namespace vcfspec::header::meta;
    other_values.reserve(other_values.size() + 4);
    other_values.emplace(struc::id, std::move(id));
    if (type == type::flag && number == number::one) number = number::zero;
    other_values.emplace(struc::number, std::move(number));
    other_values.emplace(struc::type, std::move(type));
    other_values.emplace(struc::description, add_quotes(description));
    add_structured_field(tag::info, std::move(other_values));
    return *this;
}

VcfHeader::Builder& VcfHeader::Builder::add_filter(std::string id, std::string description,
                                                   std::unordered_map<std::string, std::string> other_values)
{
    using namespace vcfspec::header::meta;
    other_values.reserve(other_values.size() + 2);
    other_values.emplace(struc::id, std::move(id));
    other_values.emplace(struc::description, add_quotes(description));
    add_structured_field(tag::filter, std::move(other_values));
    return *this;
}

VcfHeader::Builder& VcfHeader::Builder::add_format(std::string id, std::string number,
                                                   std::string type, std::string description,
                                                   std::unordered_map<std::string, std::string> other_values)
{
    using namespace vcfspec::header::meta;
    other_values.reserve(other_values.size() + 4);
    other_values.emplace(struc::id, std::move(id));
    if (type == type::flag) {
        type = type::integer; // flags not allowed in FORMAT
        if (number == number::zero) number = number::one;
    }
    other_values.emplace(struc::number, std::move(number));
    other_values.emplace(struc::type, std::move(type));
    other_values.emplace(struc::description, add_quotes(description));
    add_structured_field(tag::format, std::move(other_values));
    return *this;
}

VcfHeader::Builder& VcfHeader::Builder::add_contig(std::string id,
                                                   std::unordered_map<std::string, std::string> other_values)
{
    using namespace vcfspec::header::meta;
    other_values.emplace(struc::id, std::move(id));
    add_structured_field(tag::contig, std::move(other_values));
    return *this;
}

VcfHeader::Builder& VcfHeader::Builder::clear_info() noexcept
{
    using namespace vcfspec::header::meta;
    structured_fields_.erase(tag::info);
    return *this;
}

VcfHeader::Builder& VcfHeader::Builder::clear_format() noexcept
{
    using namespace vcfspec::header::meta;
    structured_fields_.erase(tag::format);
    samples_.clear();
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
    result.add_info("MQ", "1", "Float", "RMS mapping quality");
    result.add_info("MQ0", "1", "Integer", "Number of MAPQ == 0 reads covering this record");
    result.add_info("NS", "1", "Integer", "Number of samples with data");
    result.add_info("SB", "1", "Float", "Strand bias at this position");
    result.add_info("SOMATIC", "0", "Flag", "Indicates that the record is a somatic mutation, for cancer genomics");
    result.add_info("VALIDATED", "0", "Flag", "Validated by follow-up experiment");
    result.add_info("1000G", "0", "Flag", "Membership in 1000 Genomes");
    
    result.add_format("GT", "1", "String", "Genotype");
    result.add_format("DP", "1", "Integer", "Read depth at this position for this sample");
    result.add_format("FT", "1", "String", "Sample genotype filter indicating if this genotype was 'called'");
    result.add_format("GL", "G", "Float", "log10-scaled genotype likelihoods");
    result.add_format("GLE", "1", "Integer", "Genotype likelihoods of heterogeneous ploidy");
    result.add_format("PL", "G", "Integer", "Phred-scaled genotype likelihoods");
    result.add_format("GP", "G", "Float", "Phred-scaled genotype posterior probabilities");
    result.add_format("GQ", "1", "Integer", "Conditional genotype quality (phred-scaled)");
    result.add_format("HQ", "1", "Integer", "Haplotype qualities");
    result.add_format("PS", "1", "Integer", "Phase set");
    result.add_format("PQ", "1", "Integer", "Phasing quality");
    result.add_format("EC", "1", "Integer", "Expected alternate allele counts");
    result.add_format("MQ", "1", "Integer", "RMS mapping quality");
    result.add_format("BQ", "1", "Integer", "RMS base quality at this position");
    
    result.add_filter("PASS", "All filters passed");
    
    return result;
}

} // namespace octopus
