//
//  vcf_header.h
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

#include "vcf_type.h"

/**
 There are two types of header lines in VCF.
 
 basic:      key=value
 structured: TAG=<keyA=valueA,...,keyB=valueB>
 
 basic field 'key's must be unique. structured field 'TAG's may not be unique (e.g. INFO), but then there must
 be some unique 'key' within <>.
 */
class VcfHeader
{
public:
    VcfHeader()  = default;
    VcfHeader(std::string file_format);
    ~VcfHeader() = default;
    
    VcfHeader(const VcfHeader&)            = default;
    VcfHeader& operator=(const VcfHeader&) = default;
    VcfHeader(VcfHeader&&)                 = default;
    VcfHeader& operator=(VcfHeader&&)      = default;
    
    void set_file_format(std::string format);
    void put_sample(std::string sample);
    void put_samples(std::vector<std::string> samples);
    void put_basic_field(std::string key, std::string value); // key=value
    void put_structured_field(std::string tag, std::unordered_map<std::string, std::string> values); // TAG=<A=..,B=..>
    
    const std::string& get_file_format() const noexcept;
    unsigned get_num_samples() const noexcept;
    std::vector<std::string> get_samples() const;
    bool has_basic_field(const std::string& key) const noexcept;
    bool has_structured_field(const std::string& tag) const noexcept;
    bool has_structured_field(const std::string& tag, const std::string& key) const noexcept;
    std::vector<std::string> get_basic_field_keys() const;
    std::vector<std::string> get_structured_field_tags() const;
    const std::string& get_basic_field_value(const std::string& key) const;
    const std::string& get_structured_field_value(const std::string& tag, const std::string& id_key, const std::string& id_value, const std::string& lookup_key) const;
    
    const std::unordered_map<std::string, std::string>& get_basic_fields() const noexcept;
    std::vector<std::unordered_map<std::string, std::string>> get_structured_fields(const std::string& tag) const;
    const std::unordered_multimap<std::string, std::unordered_map<std::string, std::string>>& get_structured_fields() const noexcept;
    
    friend std::ostream& operator<<(std::ostream& os, const VcfHeader& header);
    
private:
    // required lines
    std::string file_format_;
    std::vector<std::string> samples_;
    
    // optional lines
    std::unordered_map<std::string, std::string> basic_fields_;
    std::unordered_multimap<std::string, std::unordered_map<std::string, std::string>> structured_fields_;
    //std::unordered_map<std::string, std::string> structured_field_unique_key_;
};

const std::string& get_id_field_value(const VcfHeader& header, const std::string& tag, const std::string& id_value, const std::string& lookup_key);
const std::string& get_id_field_type(const VcfHeader& header, const std::string& tag, const std::string& id_value);

VcfType get_typed_value(const VcfHeader& header, const std::string& tag, const std::string& key, const std::string& value);
VcfType get_typed_info_value(const VcfHeader& header, const std::string& key, const std::string& value);
VcfType get_typed_format_value(const VcfHeader& header, const std::string& key, const std::string& value);

std::vector<VcfType> get_typed_values(const std::string& format_key, const std::string& field_key,
                                      const std::vector<std::string>& values, const VcfHeader& header);

std::vector<VcfType> get_typed_info_values(const std::string& field_key, const std::vector<std::string>& values,
                                           const VcfHeader& header);

class VcfRecord;
unsigned get_field_cardinality(const std::string& key, const VcfRecord& record);

std::ostream& operator<<(std::ostream& os, const VcfHeader& header);

#endif /* defined(__Octopus__vcf_header__) */
