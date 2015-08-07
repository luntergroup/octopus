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
    
    unsigned get_field_cardinality(const std::string& key, const VcfRecord& record) const;
    
    VcfType get_typed_value(const std::string& format_key, const std::string& field_key, const std::string& value) const;
    // same as above but will look at all possible formats, throws if field_key is present
    // in more than one format
    VcfType get_typed_value(const std::string& field_key, const std::string& value) const;
    
    friend std::ostream& operator<<(std::ostream& os, const VcfHeader& header);
    
private:
    using Fields = std::unordered_map<std::string, std::string>;
    
    // required lines
    std::string file_format_;
    
    // optional single-lines
    Fields fields_;
    
    // optional field format lines
    std::unordered_multimap<std::string, Fields> formats_;
    
    // for fast access to types
    std::unordered_map<std::string, std::string> key_type_map_;
    
    void insert_format_line(const std::string& line);
    void insert_field_line(const std::string& line);
    void insert_lines(const std::string& lines);
};

std::ostream& operator<<(std::ostream& os, const VcfHeader& header);

#endif /* defined(__Octopus__vcf_header__) */
