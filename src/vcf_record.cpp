//
//  vcf_record.cpp
//  Octopus
//
//  Created by Daniel Cooke on 28/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "vcf_record.h"

#include <algorithm> // std::find
#include <iterator>  // std::cbegin, std::cend

const std::string& VcfRecord::get_chromosome_name() const noexcept
{
    return chrom_;
}

VcfRecord::SizeType VcfRecord::get_position() const noexcept
{
    return pos_;
}

const std::string& VcfRecord::get_id() const noexcept
{
    return id_;
}

const VcfRecord::SequenceType& VcfRecord::get_ref_allele() const noexcept
{
    return ref_;
}

unsigned VcfRecord::get_num_alt_alleles() const noexcept
{
    return static_cast<unsigned>(alt_.size());
}

const VcfRecord::SequenceType& VcfRecord::get_alt_allele(unsigned n) const noexcept
{
    return alt_.at(n);
}

VcfRecord::QualityType VcfRecord::get_quality() const noexcept
{
    return qual_;
}

bool VcfRecord::has_filter(const std::string& filter) const noexcept
{
    return std::find(std::cbegin(filter_), std::cend(filter_), filter) != std::cend(filter_);
}

unsigned VcfRecord::get_num_samples() const noexcept
{
    return static_cast<unsigned>(format_.size());
}

bool VcfRecord::has_info(const std::string& key) const noexcept
{
    return info_.count(key) == 1;
}

const std::vector<std::string>& VcfRecord::get_info_values(const std::string& key) const
{
    return info_.at(key);
}

bool VcfRecord::has_genotype_data() const noexcept
{
    return !format_.empty();
}

// non-member functions

// these INFO keys are reserved
static const std::string Info_ancestral_allele {"AA"};
static const std::string Info_genotype_allele_count {"AC"};
static const std::string Info_dbsnp {"DB"};
static const std::string Info_hapmap2 {"H2"};
static const std::string Info_hapmap3 {"H3"};
static const std::string Info_1000g {"1000G"};
static const std::string Info_somatic {"SOMATIC"};
static const std::string Info_validated {"VALIDATED"};

bool is_dbsnp_member(const VcfRecord& record) noexcept
{
    return record.has_info(Info_dbsnp);
}

bool is_hapmap2_member(const VcfRecord& record) noexcept
{
    return record.has_info(Info_hapmap2);
}

bool is_hapmap3_member(const VcfRecord& record) noexcept
{
    return record.has_info(Info_hapmap3);
}

bool is_1000g_member(const VcfRecord& record) noexcept
{
    return record.has_info(Info_1000g);
}

bool is_somatic(const VcfRecord& record) noexcept
{
    return record.has_info(Info_somatic);
}

bool is_validated(const VcfRecord& record) noexcept
{
    return record.has_info(Info_validated);
}

std::ostream& print_vector(std::ostream& os, const std::vector<std::string>& v, const std::string& delim)
{
    if (v.empty()) {
        os << ".";
    } else {
        std::copy(v.cbegin(), std::prev(v.cend()), std::ostream_iterator<std::string>(os, delim.c_str()));
        os << v.back();
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<std::string>& v)
{
    return print_vector(os, v, ",");
}

std::ostream& print_map(std::ostream& os, const std::map<std::string, std::vector<std::string>>& m, const std::string& delim)
{
    if (m.empty()) {
        os << ".";
    } else {
        auto last = std::prev(m.cend());
        std::for_each(m.cbegin(), last, [&os, delim] (const auto& p) {
            os << p.first;
            if (!p.second.empty()) {
                os << "=" << p.second;
            }
            os << delim;
        });
        os << last->first;
        if (!last->second.empty()) {
            os << "=" << last->second;
        }
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::map<std::string, std::vector<std::string>>& m)
{
    return print_map(os, m, ";");
}

std::ostream& operator<<(std::ostream& os, const VcfRecord& record)
{
    os << record.chrom_ << "\t" << record.pos_ << "\t" << record.id_ << "\t" << record.ref_ << "\t";
    os << record.alt_ << "\t";
    os << static_cast<unsigned>(record.qual_) << "\t";
    os << record.filter_ << "\t";
    os << record.info_ << "\t";
    
    if (record.has_genotype_data()) {
        print_vector(os, record.format_, ":") << "\t";
        print_map(os, record.genotypes_, ":");
    }
    
    return os;
}
