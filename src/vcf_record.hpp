//
//  vcf_record.hpp
//  Octopus
//
//  Created by Daniel Cooke on 28/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__vcf_record__
#define __Octopus__vcf_record__

#include <cstdint>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <ostream>
#include <utility>
#include <initializer_list>

#include "comparable.hpp"

// TODO: consider using #include <boost/container/small_vector.hpp> for INFO and genotype fields

class VcfRecord : public Comparable<VcfRecord>
{
public:
    class Builder;
    
    using SizeType     = std::uint_fast32_t;
    using IdType       = std::string;
    using SequenceType = std::string;
    using QualityType  = float;
    using SampleIdType = std::string;
    using KeyType      = std::string;
    using ValueType    = std::string;
    
    VcfRecord()  = default;
    
    // constructor without genotype fields
    template <typename StringType1_, typename StringType2_, typename SequenceType1_, typename SequenceType2_,
              typename Filters_, typename Info_>
    VcfRecord(StringType1_&& chrom, SizeType pos, StringType2_&& id, SequenceType1_&& ref, SequenceType2_&& alt,
              QualityType qual, Filters_&& filters, Info_&& info);
    
    // constructor with genotype fields
    template <typename StringType1_, typename StringType2_, typename SequenceType1_, typename SequenceType2_,
    typename Filters_, typename Info_, typename Format_, typename Genotypes_, typename Samples_>
    VcfRecord(StringType1_&& chrom, SizeType pos, StringType2_&& id, SequenceType1_&& ref, SequenceType2_&& alt,
              QualityType qual, Filters_&& filters, Info_&& info, Format_&& format, Genotypes_&& genotypes, Samples_&& samples);
    ~VcfRecord() = default;
    
    VcfRecord(const VcfRecord&)            = default;
    VcfRecord& operator=(const VcfRecord&) = default;
    VcfRecord(VcfRecord&&)                 = default;
    VcfRecord& operator=(VcfRecord&&)      = default;
    
    const std::string& get_chromosome_name() const noexcept;
    SizeType get_position() const noexcept;
    const IdType& get_id() const noexcept;
    const SequenceType& get_ref_allele() const noexcept;
    unsigned num_alt_alleles() const noexcept;
    const std::vector<SequenceType>& get_alt_alleles() const noexcept;
    QualityType get_quality() const noexcept;
    bool has_filter(const KeyType& filter) const noexcept;
    const std::vector<KeyType> get_filters() const noexcept;
    bool has_info(const KeyType& key) const noexcept;
    std::vector<KeyType> get_info_keys() const;
    const std::vector<ValueType>& get_info_value(const KeyType& key) const;
    
    // sample related functions
    bool has_format(const KeyType& key) const noexcept;
    unsigned format_cardinality(const KeyType& key) const noexcept;
    const std::vector<KeyType>& get_format() const noexcept;
    unsigned num_samples() const noexcept;
    bool has_genotypes() const noexcept;
    unsigned sample_ploidy() const noexcept;
    bool is_sample_phased(const SampleIdType& sample) const;
    bool is_homozygous(const SampleIdType& sample) const;
    bool is_heterozygous(const SampleIdType& sample) const;
    bool is_homozygous_ref(const SampleIdType& sample) const;
    bool is_homozygous_non_ref(const SampleIdType& sample) const;
    bool has_ref_allele(const SampleIdType& sample) const;
    bool has_alt_allele(const SampleIdType& sample) const;
    
    const std::vector<ValueType>& get_sample_value(const SampleIdType& sample, const KeyType& key) const;
    
    friend std::ostream& operator<<(std::ostream& os, const VcfRecord& record);
    
private:
    using Genotype = std::pair<std::vector<SequenceType>, bool>;
    using ValueMap = std::unordered_map<KeyType, std::vector<ValueType>>;
    
    // mandatory fields
    std::string chromosome_;
    SizeType position_;
    IdType id_;
    SequenceType ref_allele_;
    std::vector<SequenceType> alt_alleles_;
    QualityType quality_;
    std::vector<KeyType> filters_;
    ValueMap info_;
    
    // optional fields
    std::vector<KeyType> format_;
    std::unordered_map<SampleIdType, Genotype> genotypes_;
    std::unordered_map<SampleIdType, ValueMap> samples_;
    
    std::string get_allele_number(const SequenceType& allele) const;
    
    std::vector<SampleIdType> get_samples() const;
    void print_info(std::ostream& os) const;
    void print_genotype_allele_numbers(std::ostream& os, const SampleIdType& sample) const;
    void print_other_sample_data(std::ostream& os, const SampleIdType& sample) const;
    void print_sample_data(std::ostream& os) const;
};

template <typename StringType1_, typename StringType2_, typename SequenceType1_, typename SequenceType2_,
          typename Filters_, typename Info_>
VcfRecord::VcfRecord(StringType1_&& chrom, SizeType pos, StringType2_&& id, SequenceType1_&& ref,
                     SequenceType2_&& alt, QualityType qual, Filters_&& filters, Info_&& info)
:
chromosome_ {std::forward<StringType1_>(chrom)},
position_ {pos},
id_ {std::forward<StringType2_>(id)},
ref_allele_ {std::forward<SequenceType1_>(ref)},
alt_alleles_ {std::forward<SequenceType2_>(alt)},
quality_ {qual},
filters_ {std::forward<Filters_>(filters)},
info_ {std::forward<Info_>(info)},
format_ {},
genotypes_ {},
samples_ {}
{}

template <typename StringType1_, typename StringType2_, typename SequenceType1_, typename SequenceType2_,
          typename Filters_, typename Info_, typename Format_, typename Genotypes_, typename Samples_>
VcfRecord::VcfRecord(StringType1_&& chrom, SizeType pos, StringType2_&& id, SequenceType1_&& ref, SequenceType2_&& alt,
                     QualityType qual, Filters_&& filters, Info_&& info, Format_&& format, Genotypes_&& genotypes,
                     Samples_&& samples)
:
chromosome_ {std::forward<StringType1_>(chrom)},
position_ {pos},
id_ {std::forward<StringType2_>(id)},
ref_allele_ {std::forward<SequenceType1_>(ref)},
alt_alleles_ {std::forward<SequenceType2_>(alt)},
quality_ {qual},
filters_ {std::forward<Filters_>(filters)},
info_ {std::forward<Info_>(info)},
format_ {std::forward<Format_>(format)},
genotypes_ {std::forward<Genotypes_>(genotypes)},
samples_ {std::forward<Samples_>(samples)}
{}

// non-member functions

VcfRecord::SequenceType get_ancestral_allele(const VcfRecord& record);
std::vector<unsigned> get_allele_count(const VcfRecord& record);
std::vector<double> get_allele_frequency(const VcfRecord& record);
bool is_dbsnp_member(const VcfRecord& record) noexcept;
bool is_hapmap2_member(const VcfRecord& record) noexcept;
bool is_hapmap3_member(const VcfRecord& record) noexcept;
bool is_1000g_member(const VcfRecord& record) noexcept;
bool is_somatic(const VcfRecord& record) noexcept;
bool is_validated(const VcfRecord& record) noexcept;

bool operator==(const VcfRecord& lhs, const VcfRecord& rhs);
bool operator<(const VcfRecord& lhs, const VcfRecord& rhs);

namespace std {
    template <> struct hash<VcfRecord>
    {
        size_t operator()(const VcfRecord& record) const
        {
            return hash<string>()(record.get_id());
        }
    };
} // namespace std

std::ostream& operator<<(std::ostream& os, const VcfRecord& record);

class VcfRecord::Builder
{
public:
    using SizeType     = VcfRecord::SizeType;
    using IdType       = VcfRecord::IdType;
    using SequenceType = VcfRecord::SequenceType;
    using QualityType  = VcfRecord::QualityType;
    using SampleIdType = VcfRecord::SampleIdType;
    using KeyType      = VcfRecord::KeyType;
    using ValueType    = VcfRecord::ValueType;
    
    enum class Phasing { Phased, Unphased };
    
    Builder() = default;
    
    Builder& set_chromosome(std::string chromosome);
    Builder& set_position(SizeType position);
    Builder& set_id(IdType id);
    Builder& set_ref_allele(const char ref_allele);
    Builder& set_ref_allele(SequenceType ref_allele);
    Builder& set_alt_allele(const char alt_allele); // if just one
    Builder& set_alt_allele(SequenceType alt_allele); // if just one
    Builder& set_alt_alleles(std::vector<SequenceType> alt_alleles);
    Builder& set_refcall();
    Builder& set_quality(QualityType quality);
    Builder& set_filters(const std::vector<KeyType>& filters);
    Builder& set_filters(const std::initializer_list<KeyType>& filters);
    Builder& add_info(const KeyType& key, const ValueType& value);
    Builder& add_info(const KeyType& key, const std::vector<ValueType>& values);
    Builder& add_info(const KeyType& key, const std::initializer_list<ValueType>& values);
    Builder& set_format(const std::vector<KeyType>& format);
    Builder& set_format(const std::initializer_list<KeyType>& format);
    Builder& add_homozygous_ref_genotype(const SampleIdType& sample, unsigned ploidy);
    Builder& add_genotype(const SampleIdType& sample, const std::vector<SequenceType>& alleles, Phasing phasing);
    Builder& add_genotype(const SampleIdType& sample, const std::vector<unsigned>& alleles, Phasing is_phased);
    Builder& add_genotype_field(const SampleIdType& sample, const KeyType& key, const ValueType& value);
    Builder& add_genotype_field(const SampleIdType& sample, const KeyType& key, const std::vector<ValueType>& values);
    Builder& add_genotype_field(const SampleIdType& sample, const KeyType& key, const std::initializer_list<ValueType>& values);
    
    SizeType get_position() const noexcept;
    
    VcfRecord build() const;
    VcfRecord build_once() noexcept;
    
private:
    std::string chromosome_ = ".";
    SizeType position_ = 0;
    IdType id_ = ".";
    SequenceType ref_allele_ = ".";
    std::vector<SequenceType> alt_alleles_ = {"."};
    QualityType quality_ = 0;
    std::vector<KeyType> filters_ = {};
    std::unordered_map<KeyType, std::vector<ValueType>> info_ = {};
    std::vector<KeyType> format_ = {};
    std::unordered_map<SampleIdType, Genotype> genotypes_ = {};
    std::unordered_map<SampleIdType, std::unordered_map<KeyType, std::vector<ValueType>>> samples_ = {};
};

#endif /* defined(__Octopus__vcf_record__) */
