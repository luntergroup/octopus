//
//  vcf_header.hpp
//  Octopus
//
//  Created by Daniel Cooke on 28/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__vcf_header__
#define __Octopus__vcf_header__

#include <string>
#include <vector>
#include <unordered_map>
#include <ostream>

#include "equitable.hpp"
#include "vcf_type.hpp"

/*
 There are two types of header lines in VCF.
 
 basic:      key=value
 structured: TAG=<keyA=valueA,...,keyB=valueB>
 
 basic field 'key's must be unique. structured field 'TAG's may not be unique (e.g. INFO), but then there must
 be some unique 'key' within <>.
 */
class VcfHeader : public Equitable<VcfHeader>
{
public:
    class Builder;
    
    using KeyType          = std::string; // The same KeyType is used for both basic and structured fields
    using TagType          = std::string;
    using ValueType        = std::string;
    using BasicFieldMap      = std::unordered_map<KeyType, ValueType>;
    using StructuredField  = std::unordered_map<KeyType, ValueType>;
    using StructuredFieldMap = std::unordered_multimap<TagType, StructuredField>;
    
    VcfHeader()  = default;
    
    explicit VcfHeader(std::string file_format);
    template <typename T, typename U, typename B, typename S>
    explicit VcfHeader(T&& file_format, U&& samples, B&& basic_fields, S&& structured_fields);
    
    ~VcfHeader() = default;
    
    VcfHeader(const VcfHeader&)            = default;
    VcfHeader& operator=(const VcfHeader&) = default;
    VcfHeader(VcfHeader&&)                 = default;
    VcfHeader& operator=(VcfHeader&&)      = default;
    
    const ValueType& file_format() const noexcept;
    
    unsigned num_samples() const noexcept;
    std::vector<std::string> samples() const;
    
    bool has_basic_field(const KeyType& key) const noexcept;
    bool has_structured_field(const TagType& tag) const noexcept;
    bool has_structured_field(const TagType& tag, const KeyType& key) const noexcept;
    
    std::vector<KeyType> get_basic_field_keys() const;
    std::vector<TagType> get_structured_field_tags() const;
    const ValueType& get_basic_field_value(const KeyType& key) const;
    const ValueType& get_structured_field_value(const TagType& tag, const KeyType& id_key,
                                                const ValueType& id_value, const KeyType& lookup_key) const;
    
    const BasicFieldMap& get_basic_fields() const noexcept;
    std::vector<StructuredField> get_structured_fields(const TagType& tag) const;
    const StructuredFieldMap& get_structured_fields() const noexcept;
    
    friend std::ostream& operator<<(std::ostream& os, const VcfHeader& header);
    
private:
    // required fields
    std::string file_format_;
    std::vector<std::string> samples_;
    
    // optional fields
    BasicFieldMap basic_fields_;
    StructuredFieldMap structured_fields_;
};

template <typename T, typename U, typename B, typename S>
VcfHeader::VcfHeader(T&& file_format, U&& samples, B&& basic_fields, S&& structured_fields)
:
file_format_ {std::forward<T>(file_format)},
samples_ {std::forward<U>(samples)},
basic_fields_ {std::forward<B>(basic_fields)},
structured_fields_ {std::forward<S>(structured_fields)}
{}

const std::string& get_id_field_value(const VcfHeader& header, const VcfHeader::TagType& tag,
                                      const VcfHeader::ValueType& id_value, const VcfHeader::KeyType& lookup_key);
const std::string& get_id_field_type(const VcfHeader& header, const VcfHeader::TagType& tag,
                                     const VcfHeader::ValueType& id_value);

VcfType get_typed_value(const VcfHeader& header, const VcfHeader::TagType& tag,
                        const VcfHeader::KeyType& key, const VcfHeader::ValueType& value);

VcfType get_typed_info_value(const VcfHeader& header, const VcfHeader::KeyType& key,
                             const VcfHeader::ValueType& value);

VcfType get_typed_format_value(const VcfHeader& header, const VcfHeader::KeyType& key,
                               const VcfHeader::ValueType& value);

std::vector<VcfType> get_typed_values(const VcfHeader& header, const VcfHeader::KeyType& format_key,
                                      const VcfHeader::KeyType& field_key,
                                      const std::vector<VcfHeader::ValueType>& values);

std::vector<VcfType> get_typed_info_values(const VcfHeader& header, const VcfHeader::KeyType& field_key,
                                           const std::vector<VcfHeader::ValueType>& values);

std::vector<VcfType> get_typed_format_values(const VcfHeader& header, const VcfHeader::KeyType& field_key,
                                             const std::vector<VcfHeader::ValueType>& values);

bool contig_line_exists(const VcfHeader& header, const std::string& contig);

bool operator==(const VcfHeader& lhs, const VcfHeader& rhs);

namespace std {
    template <> struct hash<VcfHeader>
    {
        size_t operator()(const VcfHeader& header) const
        {
            return hash<VcfHeader::ValueType>()(header.file_format());
        }
    };
} // namespace std

std::ostream& operator<<(std::ostream& os, const VcfHeader& header);

class VcfHeader::Builder
{
public:
    Builder() = default;
    Builder(const VcfHeader& header);
    
    Builder& set_file_format(std::string file_format);
    Builder& add_sample(std::string sample);
    Builder& set_samples(std::vector<std::string> samples);
    Builder& add_basic_field(std::string key, std::string value);
    Builder& add_structured_field(std::string tag, std::unordered_map<std::string, std::string> values);
    
    Builder& add_info(std::string id, std::string number, std::string type, std::string description,
                      std::unordered_map<std::string, std::string> other_values = {});
    Builder& add_filter(std::string id, std::string description,
                        std::unordered_map<std::string, std::string> other_values = {});
    Builder& add_format(std::string id, std::string number, std::string type, std::string description,
                        std::unordered_map<std::string, std::string> other_values = {});
    Builder& add_contig(std::string id, std::unordered_map<std::string, std::string> other_values = {});
    
    VcfHeader build() const;
    VcfHeader build_once() noexcept;
    
private:
    std::string file_format_ = "VCFv4.2";
    std::vector<std::string> samples_ = {};
    std::unordered_map<std::string, std::string> basic_fields_ = {};
    std::unordered_multimap<std::string, std::unordered_map<std::string, std::string>> structured_fields_ = {};
};

// A VcfHeader::Builder pre-filled with all reserved INFO and FORMAT fields
VcfHeader::Builder get_default_header_builder();

#endif /* defined(__Octopus__vcf_header__) */
