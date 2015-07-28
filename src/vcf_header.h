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
#include <unordered_map>
#include <boost/optional.hpp>

using boost::optional;

class BaseHeaderLine {};

class IdHeaderLine : public BaseHeaderLine
{
public:
    const std::string& get_id() const noexcept { return id_; }
private:
    std::string id_;
};

class DescriptiveHeaderLine : public IdHeaderLine
{
public:
    const std::string& get_description() const noexcept;
private:
    std::string description_;
};

enum class HeaderLineType { Integer, Float, Flag, Chracter, String };

class BaseHeaderTypedLine : public BaseHeaderLine
{
public:
    virtual ~BaseHeaderTypedLine() = default;
    virtual unsigned get_number() const noexcept = 0;
private:
    HeaderLineType type_;
};

class InfoFormatLine : public BaseHeaderTypedLine, public DescriptiveHeaderLine
{
public:
    virtual ~InfoFormatLine() = default;
    virtual unsigned get_number() const noexcept = 0;
private:
    optional<std::string> source_;
    optional<std::string> version_;
};

class NumericInfoFormatLine : public InfoFormatLine
{
public:
    virtual unsigned get_number() const noexcept { return number_; }
private:
    unsigned number_;
};

class FlaggedInfoFormatLine : public InfoFormatLine
{
public:
    virtual unsigned get_number() const noexcept;
private:
    enum class NumberFlags { A, R, G, N };
    NumberFlags number_flag_;
};

class FilterFormatLine : public DescriptiveHeaderLine {};

class FortmatFormatLine
{
    
};

class AltAlleleFormatLine : public DescriptiveHeaderLine
{
    
};

class ContigFormatLine : public IdHeaderLine
{
public:
    const std::string get_url() const noexcept;
private:
    std::string url_;
};

class SampleFormatLine {};

class PedigreeFormatLine {};

class VcfHeader
{
public:
    VcfHeader()  = default;
    ~VcfHeader() = default;
    
    VcfHeader(const VcfHeader&)            = default;
    VcfHeader& operator=(const VcfHeader&) = default;
    VcfHeader(VcfHeader&&)                 = default;
    VcfHeader& operator=(VcfHeader&&)      = default;
    
private:
    // required lines
    std::string file_format_;
    
    // optional field format lines
    std::vector<InfoFormatLine> info_format_lines_;
    std::vector<FilterFormatLine> filter_format_lines_;
    std::vector<FortmatFormatLine> format_format_lines_;
    std::vector<AltAlleleFormatLine> alt_allele_format_lines_;
    std::vector<ContigFormatLine> contig_format_lines_;
    std::vector<SampleFormatLine> sample_format_lines_;
    std::vector<PedigreeFormatLine> pedigree_format_lines_;
    
    // optional single-lines
    optional<std::string> assembly_;
    optional<std::tm> file_date_;
    optional<std::string> source_;
    optional<std::string> reference_;
    optional<std::string> contig_;
    optional<std::string> phasing_;
    
    
};

#endif /* defined(__Octopus__vcf_header__) */
