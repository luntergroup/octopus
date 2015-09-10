//
//  vcf_header.cpp
//  Octopus
//
//  Created by Daniel Cooke on 28/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "vcf_header.h"

#include <algorithm> // std::find_if, std::transform, std::sort, std::unique
#include <utility>   // std::move
#include <iterator>  // std::back_inserter

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
    
    std::sort(result.begin(), result.end());
    
    result.erase(std::unique(result.begin(), result.end()), result.end());
    
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
    auto er = structured_fields_.equal_range(tag);
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
    
    auto er = structured_fields_.equal_range(tag);
    std::transform(er.first, er.second, std::back_inserter(result), [] (const auto& p) { return p.second; });
    
    return result;
}

const std::unordered_multimap<VcfHeader::TagType, VcfHeader::StructuredField>& VcfHeader::get_structured_fields() const noexcept
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

unsigned get_field_cardinality(const VcfHeader::KeyType& key, const VcfRecord& record)
{
    return 0;
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

VcfHeader::Builder& VcfHeader::Builder::add_structured_field(std::string tag, std::unordered_map<std::string, std::string> values)
{
    structured_fields_.emplace(std::move(tag), std::move(values));
    return *this;
}

VcfHeader::Builder& VcfHeader::Builder::add_info(std::string id, std::string number, std::string type, std::string description,
                                                 std::unordered_map<std::string, std::string> other_values)
{
    other_values.emplace("ID", std::move(id));
    other_values.emplace("Number", std::move(number));
    other_values.emplace("Type", std::move(type));
    other_values.emplace("Description", std::move(description));
    
    structured_fields_.emplace("INFO", std::move(other_values));
    
    return *this;
}

VcfHeader::Builder& VcfHeader::Builder::add_filter(std::string id, std::string description,
                                                   std::unordered_map<std::string, std::string> other_values)
{
    other_values.emplace("ID", std::move(id));
    other_values.emplace("Description", std::move(description));
    
    structured_fields_.emplace("FILTER", std::move(other_values));
    
    return *this;
}

VcfHeader::Builder& VcfHeader::Builder::add_format(std::string id, std::string number, std::string type, std::string description,
                                                   std::unordered_map<std::string, std::string> other_values)
{
    other_values.emplace("ID", std::move(id));
    other_values.emplace("Number", std::move(number));
    other_values.emplace("Type", std::move(type));
    other_values.emplace("Description", std::move(description));
    
    structured_fields_.emplace("FORMAT", std::move(other_values));
    
    return *this;
}

VcfHeader::Builder& VcfHeader::Builder::add_contig(std::string id, std::unordered_map<std::string, std::string> other_values)
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
