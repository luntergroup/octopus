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
 
 1) key=value
 2) TAG=<keyA=valueA,...,keyB=valueB>
 
 for (1), 'key' must be unique. for (2), 'TAG' may not be unique (e.g. INFO), but then there must
 be some unique 'key' within <>.
 */

class VcfRecord;

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
    void put_field(std::string key, std::string value); // key=value
    void put_field(std::string tag, std::unordered_map<std::string, std::string> values); // TAG=<A=..,B=..>
    
    const std::string& get_file_format() const noexcept;
    unsigned get_num_samples() const noexcept;
    std::vector<std::string> get_samples() const;
    bool has_field(const std::string& key) const noexcept;
    bool has_tag(const std::string& tag) const noexcept;
    bool has_field(const std::string& tag, const std::string& key) const noexcept;
    const std::string& get_field(const std::string& key) const;
    const std::string& get_field(const std::string& tag, const std::string& id_key, const std::string& id_value, const std::string& lookup_key) const;
    
    friend std::ostream& operator<<(std::ostream& os, const VcfHeader& header);
    
private:
    using Fields = std::unordered_map<std::string, std::string>;
    
    // required lines
    std::string file_format_;
    std::vector<std::string> samples_;
    
    // optional lines
    Fields fields_;
    std::unordered_multimap<std::string, Fields> formats_;
};

const std::string& get_id_field(const VcfHeader& header, const std::string& tag, const std::string& id_value, const std::string& lookup_key);
const std::string& get_id_field_type(const VcfHeader& header, const std::string& tag, const std::string& id_value);

VcfType get_typed_value(const VcfHeader& header, const std::string& tag, const std::string& key, const std::string& value);
VcfType get_typed_info_value(const VcfHeader& header, const std::string& key, const std::string& value);
VcfType get_typed_format_value(const VcfHeader& header, const std::string& key, const std::string& value);

std::vector<VcfType> get_typed_values(const std::string& format_key, const std::string& field_key,
                                      const std::vector<std::string>& values, const VcfHeader& header);

std::vector<VcfType> get_typed_info_values(const std::string& field_key, const std::vector<std::string>& values,
                                           const VcfHeader& header);

unsigned get_field_cardinality(const std::string& key, const VcfRecord& record);

std::ostream& operator<<(std::ostream& os, const VcfHeader& header);

#endif /* defined(__Octopus__vcf_header__) */
