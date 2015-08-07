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
#include <ctime>
#include <ostream>
#include <unordered_map>
#include <stdexcept>
#include <boost/optional.hpp>

#include "vcf_type.h"

using boost::optional;

class VcfRecord;

class VcfHeaderLine {};

enum class VcfHeaderLineType {Integer, Float, Character, String, Flag};

class TypeField
{
public:
    
private:
    
};

class VcfHeader
{
public:
    VcfHeader()  = default;
    VcfHeader(const std::string& lines);
    ~VcfHeader() = default;
    
    VcfHeader(const VcfHeader&)            = default;
    VcfHeader& operator=(const VcfHeader&) = default;
    VcfHeader(VcfHeader&&)                 = default;
    VcfHeader& operator=(VcfHeader&&)      = default;
    
    void insert_line(const std::string& line);
    
    const std::string& get_file_format() const noexcept;
    unsigned get_num_samples() const noexcept;
    std::vector<std::string> get_samples() const;
    
    // acessors for single fields
    const std::string& get_field(const std::string& field_key) const;
    
    unsigned get_field_cardinality(const std::string& key, const VcfRecord& record) const;
    
    // The second version will look at all format fields and throw if field_key is not unique
    VcfType get_typed_value(const std::string& format_key, const std::string& field_key, const std::string& value) const;
    VcfType get_typed_value(const std::string& field_key, const std::string& value) const;
    
    friend std::ostream& operator<<(std::ostream& os, const VcfHeader& header);
    
private:
    using Fields = std::unordered_map<std::string, std::string>;
    
    // required lines
    std::string file_format_;
    std::vector<std::string> samples_;
    
    // optional single-lines
    Fields fields_;
    
    // optional field format lines
    std::unordered_multimap<std::string, Fields> formats_;
    
    void insert_format_line(const std::string& line);
    void insert_field_line(const std::string& line);
    void insert_lines(const std::string& lines);
};

std::ostream& operator<<(std::ostream& os, const VcfHeader& header);

#endif /* defined(__Octopus__vcf_header__) */
