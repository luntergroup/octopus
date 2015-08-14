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

unsigned VcfRecord::num_alt_alleles() const noexcept
{
    return static_cast<unsigned>(alt_alleles_.size());
}

const std::vector<VcfRecord::SequenceType>& VcfRecord::get_alt_alleles() const noexcept
{
    return alt_alleles_;
}

VcfRecord::QualityType VcfRecord::get_quality() const noexcept
{
    return quality_;
}

bool VcfRecord::has_filter(const KeyType& filter) const noexcept
{
    return std::find(std::cbegin(filters_), std::cend(filters_), filter) != std::cend(filters_);
}

const std::vector<VcfRecord::KeyType> VcfRecord::get_filters() const noexcept
{
    return filters_;
}

bool VcfRecord::has_info(const KeyType& key) const noexcept
{
    return info_.count(key) == 1;
}

std::vector<VcfRecord::KeyType> VcfRecord::get_info_keys() const
{
    std::vector<KeyType> result {};
    result.reserve(info_.size());
    
    std::transform(info_.cbegin(), info_.cend(), std::back_inserter(result), [] (const auto& p) {
        return p.first;
    });
    
    return result;
}

const std::vector<std::string>& VcfRecord::get_info_value(const KeyType& key) const
{
    return info_.at(key);
}

bool VcfRecord::has_format(const KeyType& key) const noexcept
{
    return std::find(std::cbegin(format_), std::cend(format_), key) != std::cend(format_);
}

unsigned VcfRecord::format_cardinality(const KeyType& key) const noexcept
{
    return (has_format(key)) ? static_cast<unsigned>(samples_.cbegin()->second.at(key).size()) : 0;
}

const std::vector<VcfRecord::KeyType>& VcfRecord::get_format() const noexcept
{
    return format_;
}

unsigned VcfRecord::num_samples() const noexcept
{
    return static_cast<unsigned>((has_genotypes()) ? genotypes_.size() : samples_.size());
}

bool VcfRecord::has_genotypes() const noexcept
{
    return !genotypes_.empty();
}

unsigned VcfRecord::sample_ploidy() const noexcept
{
    // all samples must have the same ploidy
    return (has_genotypes()) ? static_cast<unsigned>(genotypes_.cbegin()->second.first.size()) : 0;
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

const std::vector<std::string>& VcfRecord::get_sample_value(const SampleIdType& sample, const KeyType& key) const
{
    return (key == "GT") ? genotypes_.at(sample).first : samples_.at(sample).at(key);
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
    if (num_samples() > 0) {
        //print_vector(os, format_, ":");
        os << "\t";
        
        bool has_genotype {has_genotypes()};
        
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
    os << record.chromosome_ << "\t";
    os << record.position_ << "\t";
    os << record.id_ << "\t";
    os << record.ref_allele_ << "\t";
    os << record.alt_alleles_ << "\t";
    os << static_cast<unsigned>(record.quality_) << "\t";
    os << record.filters_ << "\t";
    record.print_info(os);
    os << "\t";
    record.print_sample_data(os);
    
    return os;
}

// VcfRecord::Builder

VcfRecord::Builder& VcfRecord::Builder::set_chromosome(const std::string& chromosome)
{
    chromosome_ = chromosome;
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_position(SizeType position)
{
    position_ = position;
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_id(const std::string& id)
{
    id_ = id;
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_ref_allele(const SequenceType& ref_allele)
{
    ref_allele_ = ref_allele;
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_alt_allele(const SequenceType& alt_allele)
{
    alt_alleles_[0] = alt_allele;
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_alt_alleles(const std::vector<SequenceType>& alt_alleles)
{
    alt_alleles_ = alt_alleles;
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_quality(QualityType quality)
{
    quality_ = quality;
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_filters(const std::vector<KeyType>& filters)
{
    filters_ = filters;
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::add_info(const KeyType& key, const std::vector<std::string>& values)
{
    info_.emplace(key, values);
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_format(const std::vector<KeyType>& format)
{
    format_ = format;
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::add_genotype(const SampleIdType& sample, const std::vector<SequenceType>& alleles, bool is_phased)
{
    genotypes_.emplace(sample, std::make_pair(alleles, is_phased));
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::add_genotype_field(const SampleIdType& sample, const KeyType& key, const std::vector<std::string>& values)
{
    samples_[sample].emplace(key, values);
    return *this;
}

VcfRecord VcfRecord::Builder::build() const
{
    if (genotypes_.empty() && samples_.empty()) {
        return VcfRecord {chromosome_, position_, id_, ref_allele_, alt_alleles_, quality_, filters_, info_};
    } else {
        return VcfRecord {chromosome_, position_, id_, ref_allele_, alt_alleles_, quality_, filters_, info_, format_, genotypes_, samples_};
    }
}

VcfRecord VcfRecord::Builder::build_once() noexcept
{
    if (genotypes_.empty() && samples_.empty()) {
        return VcfRecord {std::move(chromosome_), position_, std::move(id_), std::move(ref_allele_),
            std::move(alt_alleles_), quality_, std::move(filters_), std::move(info_)};
    } else {
        return VcfRecord {std::move(chromosome_), position_, std::move(id_), std::move(ref_allele_),
            std::move(alt_alleles_), quality_, std::move(filters_), std::move(info_),
            std::move(format_), std::move(genotypes_), std::move(samples_)};
    }
}
