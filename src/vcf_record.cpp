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
    return chromosome_;
}

VcfRecord::SizeType VcfRecord::get_position() const noexcept
{
    return position_;
}

const std::string& VcfRecord::get_id() const noexcept
{
    return id_;
}

const VcfRecord::SequenceType& VcfRecord::get_ref_allele() const noexcept
{
    return ref_allele_;
}

unsigned VcfRecord::get_num_alt_alleles() const noexcept
{
    return static_cast<unsigned>(alt_alleles_.size());
}

const VcfRecord::SequenceType& VcfRecord::get_alt_allele(unsigned n) const noexcept
{
    return alt_alleles_.at(n);
}

VcfRecord::QualityType VcfRecord::get_quality() const noexcept
{
    return quality_;
}

bool VcfRecord::has_filter(const KeyType& filter) const noexcept
{
    return std::find(std::cbegin(filters_), std::cend(filters_), filter) != std::cend(filters_);
}

bool VcfRecord::has_info(const KeyType& key) const noexcept
{
    return info_.count(key) == 1;
}

const std::vector<std::string>& VcfRecord::get_info_value(const KeyType& key) const
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

bool VcfRecord::is_sample_phased(const SampleIdType& sample) const
{
    return genotypes_.at(sample).second;
}

bool VcfRecord::is_homozygous(const SampleIdType& sample) const
{
    const auto& genotype = genotypes_.at(sample).first;
    return std::adjacent_find(std::cbegin(genotype), std::cend(genotype),
                              std::not_equal_to<SequenceType>()) == std::cend(genotype);
}

bool VcfRecord::is_heterozygous(const SampleIdType& sample) const
{
    return !is_homozygous(sample);
}

bool VcfRecord::is_homozygous_ref(const SampleIdType& sample) const
{
    const auto& genotype = genotypes_.at(sample).first;
    return std::all_of(std::cbegin(genotype), std::cend(genotype),
                       [this] (const auto& allele) { return allele == ref_allele_; });
}

bool VcfRecord::is_homozygous_non_ref(const SampleIdType& sample) const
{
    const auto& genotype = genotypes_.at(sample).first;
    return genotype.front() != ref_allele_ && is_homozygous(sample);
}

bool VcfRecord::has_ref_allele(const SampleIdType& sample) const
{
    const auto& genotype = genotypes_.at(sample).first;
    return std::find(std::cbegin(genotype), std::cend(genotype), ref_allele_) != std::cend(genotype);
}

bool VcfRecord::has_alt_allele(const SampleIdType& sample) const
{
    const auto& genotype = genotypes_.at(sample).first;
    return std::find_if_not(std::cbegin(genotype), std::cend(genotype),
                            [this] (const auto& allele) { return allele == ref_allele_; }) != std::cend(genotype);
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

std::string VcfRecord::get_allele_number(const SequenceType& allele) const
{
    if (allele == ".") {
        return ".";
    } else if (allele == ref_allele_) {
        return "0";
    } else {
        auto it = std::find(std::cbegin(alt_alleles_), std::cend(alt_alleles_), allele);
        return std::to_string(std::distance(std::cbegin(alt_alleles_), it) + 1);
    }
}

void VcfRecord::print_info(std::ostream& os) const
{
    if (info_.empty()) {
        os << ".";
    } else {
        auto last = std::next(info_.cbegin(), info_.size() - 1);
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
                   [this] (const auto& allele) { return get_allele_number(allele); });
    print_vector(os, allele_numbers, (genotype.second) ? "|" : "/");
}

void VcfRecord::print_other_sample_data(std::ostream& os, const SampleIdType& sample) const
{
    const auto& data = samples_.at(sample);
    if (data.empty()) {
        os << ".";
    } else {
        auto last = std::next(data.cbegin(), data.size() - 1);
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
        //print_vector(os, format_, ":");
        os << "\t";
        
        bool has_genotype {has_genotype_data()};
        
        auto last = std::next(samples_.cbegin(), samples_.size() - 1);
        std::for_each(std::cbegin(samples_), last, [this, &os, has_genotype] (const auto& sample_data) {
            if (has_genotype) {
                print_genotype_allele_numbers(os, sample_data.first);
                if (!samples_.empty()) {
                    os << ":";
                }
            }
            print_other_sample_data(os, sample_data.first);
            os << "\t";
        });
        
        const auto& sample_data = *last;
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
static const VcfRecord::KeyType Info_ancestral_allele {"AA"};
static const VcfRecord::KeyType Info_genotype_allele_count {"AC"};
static const VcfRecord::KeyType Info_dbsnp {"DB"};
static const VcfRecord::KeyType Info_hapmap2 {"H2"};
static const VcfRecord::KeyType Info_hapmap3 {"H3"};
static const VcfRecord::KeyType Info_1000g {"1000G"};
static const VcfRecord::KeyType Info_somatic {"SOMATIC"};
static const VcfRecord::KeyType Info_validated {"VALIDATED"};

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
//    os << record.chromosome_ << "\t";
//    os << record.position_ << "\t";
//    os << record.id_ << "\t";
//    os << record.ref_allele_ << "\t";
//    os << record.alt_alleles_ << "\t";
//    os << static_cast<unsigned>(record.quality_) << "\t";
//    os << record.filters_ << "\t";
//    record.print_info(os);
//    os << "\t";
//    record.print_sample_data(os);
    
    return os;
}
