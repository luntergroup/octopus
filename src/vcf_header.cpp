//
//  vcf_header.cpp
//  Octopus
//
//  Created by Daniel Cooke on 28/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "vcf_header.h"

#include <algorithm> // std::find_if, std::transform
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

void VcfHeader::set_file_format(std::string format)
{
    file_format_ = std::move(format);
}

void VcfHeader::put_sample(std::string sample)
{
    samples_.emplace_back(std::move(sample));
}

void VcfHeader::put_samples(std::vector<std::string> samples)
{
    samples_.insert(samples_.end(), std::make_move_iterator(samples.begin()), std::make_move_iterator(samples.end()));
}

void VcfHeader::put_field(std::string key, std::string value)
{
    if (key != "fileformat") {
        fields_.emplace(std::move(key), std::move(value));
    }
}

void VcfHeader::put_field(std::string tag, std::unordered_map<std::string, std::string> values)
{
    formats_.emplace(std::move(tag), std::move(values));
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

bool VcfHeader::has_field(const std::string& key) const noexcept
{
    return fields_.count(key) == 1;
}

bool VcfHeader::has_tag(const std::string& tag) const noexcept
{
    return formats_.count(tag) > 0;
}

bool VcfHeader::has_field(const std::string& tag, const std::string& key) const noexcept
{
    return formats_.count(tag) > 0; // TODO: complete this
}

const std::string& VcfHeader::get_field(const std::string& key) const
{
    return fields_.at(key);
}

const std::string& VcfHeader::get_field(const std::string& tag, const std::string& id_key,
                                        const std::string& id_value, const std::string& lookup_key) const
{
    auto er = formats_.equal_range(tag);
    return std::find_if(er.first, er.second,
                        [&id_key, &id_value, &lookup_key] (const auto& p) {
                            return p.second.at(id_key) == id_value;
                        })->second.at(lookup_key);
}

// private methods

// non-member methods

const std::string& get_id_field(const VcfHeader& header, const std::string& tag, const std::string& id_value, const std::string& lookup_key)
{
    return header.get_field(tag, "ID", id_value, lookup_key);
}

const std::string& get_id_field_type(const VcfHeader& header, const std::string& tag, const std::string& id_value)
{
    return header.get_field(tag, "ID", id_value, "Type");
}

VcfType get_typed_value(const VcfHeader& header, const std::string& tag, const std::string& key, const std::string& value)
{
    return make_vcf_type(get_id_field_type(header, tag, key), value);
}

VcfType get_typed_info_value(const VcfHeader& header, const std::string& key, const std::string& value)
{
    return get_typed_value(header, "INFO", key, value);
}

VcfType get_typed_format_value(const VcfHeader& header, const std::string& key, const std::string& value)
{
    return get_typed_value(header, "FORMAT", key, value);
}

std::vector<VcfType> get_typed_values(const std::string& format_key, const std::string& field_key,
                                      const std::vector<std::string>& values, const VcfHeader& header)
{
    std::vector<VcfType> result {};
    result.reserve(values.size());
    std::transform(values.cbegin(), values.cend(), std::back_inserter(result),
                   [&header, &format_key, &field_key] (const auto& value) {
                       return get_typed_value(header, format_key, field_key, value);
                   });
    return result;
}

std::vector<VcfType> get_typed_info_values(const std::string& field_key, const std::vector<std::string>& values,
                                           const VcfHeader& header)
{
    return get_typed_values("INFO", field_key, values, header);
}

unsigned get_field_cardinality(const std::string& key, const VcfRecord& record)
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
    
    for (const auto& field : header.fields_) {
        os << "##" << field.first << "=" << field.second << std::endl;
    }
    
    for (const auto& format : header.formats_) {
        os << "##" << format.first << "=" << format.second << std::endl;
    }
    
    return os;
}
