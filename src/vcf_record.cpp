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

// public methods

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

bool VcfRecord::has_info(const std::string& key) const noexcept
{
    return info_.count(key) == 1;
}

const std::vector<std::string>& VcfRecord::get_info_values(const std::string& key) const
{
    return info_.at(key);
}

bool VcfRecord::has_sample_data() const noexcept
{
    return !format_.empty();
}

unsigned VcfRecord::num_samples() const noexcept
{
    return static_cast<unsigned>((has_genotype_data() ? genotypes_.size() : samples_.size()));
}

bool VcfRecord::has_genotype_data() const noexcept
{
    return !genotypes_.empty();
}

unsigned VcfRecord::sample_ploidy() const noexcept
{
    // all samples must have the same ploidy
    return (has_genotype_data()) ? static_cast<unsigned>(genotypes_.cbegin()->second.first.size()) : 0;
}

// helper non-members needed for printing

void print_vector(std::ostream& os, const std::vector<std::string>& v,
                  const std::string& delim = ",", const std::string& empty_value = ".")
{
    if (v.empty()) {
        os << empty_value;
    } else {
        std::copy(v.cbegin(), std::prev(v.cend()), std::ostream_iterator<std::string>(os, delim.c_str()));
        os << v.back();
    }
}

std::ostream& operator<<(std::ostream& os, const std::vector<std::string>& v)
{
    print_vector(os, v);
    return os;
}

// private methods

std::string VcfRecord::to_number(const SequenceType& allele) const
{
    if (allele == ".") {
        return allele;
    } else if (allele == ref_) {
        return "0";
    } else {
        return std::to_string(std::distance(std::cbegin(alt_), std::find(std::cbegin(alt_), std::cend(alt_), allele)) + 1);
    }
}

void VcfRecord::print_info(std::ostream& os) const
{
    if (info_.empty()) {
        os << ".";
    } else {
        auto last = std::prev(info_.cend());
        std::for_each(info_.cbegin(), last, [&os] (const auto& p) {
            os << p.first;
            if (!p.second.empty()) {
                os << "=" << p.second;
            }
            os << ';';
        });
        os << last->first;
        if (!last->second.empty()) {
            os << "=" << last->second;
        }
    }
}

void VcfRecord::print_genotype_allele_numbers(std::ostream& os, const SampleIdType& sample) const
{
    std::vector<std::string> allele_numbers(sample_ploidy());
    const auto& genotype = genotypes_.at(sample);
    std::transform(std::cbegin(genotype.first), std::cend(genotype.first), std::begin(allele_numbers),
                   [this] (const auto& allele) { return to_number(allele); });
    print_vector(os, allele_numbers, (genotype.second) ? "|" : "/");
}

void VcfRecord::print_other_sample_data(std::ostream& os, const SampleIdType& sample) const
{
    const auto& data = samples_.at(sample);
    if (data.empty()) {
        os << ".";
    } else {
        auto last = std::prev(data.cend());
        std::for_each(data.cbegin(), last, [&os] (const auto& p) {
            print_vector(os, p.second, ",");
            os << ":";
        });
        print_vector(os, last->second, ",");
    }
}

void VcfRecord::print_sample_data(std::ostream& os) const
{
    if (has_sample_data()) {
        print_vector(os, format_, ":");
        os << "\t";
        
        bool has_genotype {has_genotype_data()};
        
        std::for_each(std::cbegin(samples_), std::prev(std::cend(samples_)), [this, &os, has_genotype] (const auto& sample_data) {
            if (has_genotype) {
                print_genotype_allele_numbers(os, sample_data.first);
                if (!samples_.empty()) {
                    os << ":";
                }
            }
            print_other_sample_data(os, sample_data.first);
            os << "\t";
        });
        
        const auto& sample_data = *std::prev(std::cend(samples_));
        if (has_genotype) {
            print_genotype_allele_numbers(os, sample_data.first);
            if (!samples_.empty()) {
                os << ":";
            }
        }
        print_other_sample_data(os, sample_data.first);
    }
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

std::ostream& operator<<(std::ostream& os, const VcfRecord& record)
{
    os << record.chrom_ << "\t" << record.pos_ << "\t" << record.id_ << "\t" << record.ref_ << "\t";
    os << record.alt_ << "\t";
    os << static_cast<unsigned>(record.qual_) << "\t";
    os << record.filter_ << "\t";
    record.print_info(os);
    os << "\t";
    record.print_sample_data(os);
    
    return os;
}
