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
#include <functional>

#include <boost/optional.hpp>
#include <boost/container/flat_map.hpp>

#include <concepts/comparable.hpp>

namespace octopus {

// TODO: consider using #include <boost/container/small_vector.hpp> for INFO and genotype fields

class VcfRecord : public Comparable<VcfRecord>
{
public:
    class Builder;
    
    using SizeType           = std::uint_fast32_t;
    using IdType             = std::string;
    using NucleotideSequence = std::string;
    using QualityType        = float;
    using SampleName         = std::string;
    using KeyType            = std::string;
    using ValueType          = std::string;
    
    VcfRecord()  = default;
    
    // constructor without genotype fields
    template <typename String1, typename String2, typename Sequence1, typename Sequence2,
              typename Filters, typename Info>
    VcfRecord(String1&& chrom, SizeType pos, String2&& id, Sequence1&& ref, Sequence2&& alt,
              boost::optional<QualityType> qual, Filters&& filters, Info&& info);
    
    // constructor with genotype fields
    template <typename String1, typename String2, typename Sequence1, typename Sequence2,
    typename Filters, typename Info, typename Format, typename Genotypes, typename Samples>
    VcfRecord(String1&& chrom, SizeType pos, String2&& id, Sequence1&& ref, Sequence2&& alt,
              boost::optional<QualityType> qual, Filters&& filters, Info&& info,
              Format&& format, Genotypes&& genotypes, Samples&& samples);
    
    VcfRecord(const VcfRecord&)            = default;
    VcfRecord& operator=(const VcfRecord&) = default;
    VcfRecord(VcfRecord&&)                 = default;
    VcfRecord& operator=(VcfRecord&&)      = default;
    
    ~VcfRecord() = default;
    
    const std::string& chrom() const noexcept;
    SizeType pos() const noexcept;
    const IdType& id() const noexcept;
    const NucleotideSequence& ref() const noexcept;
    unsigned num_alt() const noexcept;
    const std::vector<NucleotideSequence>& alt() const noexcept;
    boost::optional<QualityType> qual() const noexcept;
    bool has_filter(const KeyType& filter) const noexcept;
    const std::vector<KeyType> filter() const noexcept;
    bool has_info(const KeyType& key) const noexcept;
    std::vector<KeyType> info_keys() const;
    const std::vector<ValueType>& info_value(const KeyType& key) const;
    
    // sample related functions
    bool has_format(const KeyType& key) const noexcept;
    unsigned format_cardinality(const KeyType& key) const noexcept;
    const std::vector<KeyType>& format() const noexcept;
    unsigned num_samples() const noexcept;
    bool has_genotypes() const noexcept;
    unsigned ploidy(const SampleName& sample) const;
    bool is_sample_phased(const SampleName& sample) const;
    bool is_homozygous(const SampleName& sample) const;
    bool is_heterozygous(const SampleName& sample) const;
    bool is_homozygous_ref(const SampleName& sample) const;
    bool is_homozygous_non_ref(const SampleName& sample) const;
    bool has_ref_allele(const SampleName& sample) const;
    bool has_alt_allele(const SampleName& sample) const;
    
    const std::vector<ValueType>& get_sample_value(const SampleName& sample, const KeyType& key) const;
    
    friend std::ostream& operator<<(std::ostream& os, const VcfRecord& record);
    
    friend Builder;
    
private:
    using Genotype = std::pair<std::vector<NucleotideSequence>, bool>;
    using ValueMap = boost::container::flat_map<KeyType, std::vector<ValueType>>;
    
    // mandatory fields
    std::string chrom_;
    SizeType pos_;
    IdType id_;
    NucleotideSequence ref_;
    std::vector<NucleotideSequence> alt_;
    boost::optional<QualityType> qual_;
    std::vector<KeyType> filter_;
    ValueMap info_;
    
    // optional fields
    std::vector<KeyType> format_;
    boost::container::flat_map<SampleName, Genotype> genotypes_;
    boost::container::flat_map<SampleName, ValueMap> samples_;
    
    std::string get_allele_number(const NucleotideSequence& allele) const;
    
    std::vector<SampleName> samples() const;
    void print_info(std::ostream& os) const;
    void print_genotype_allele_numbers(std::ostream& os, const SampleName& sample) const;
    void print_other_sample_data(std::ostream& os, const SampleName& sample) const;
    void print_sample_data(std::ostream& os) const;
};

template <typename String1, typename String2, typename Sequence1, typename Sequence2,
          typename Filters, typename Info>
VcfRecord::VcfRecord(String1&& chrom, SizeType pos, String2&& id,  Sequence1&& ref, Sequence2&& alt,
                     boost::optional<QualityType> qual, Filters&& filters, Info&& info)
:
chrom_ {std::forward<String1>(chrom)},
pos_ {pos},
id_ {std::forward<String2>(id)},
ref_ {std::forward<Sequence1>(ref)},
alt_ {std::forward<Sequence2>(alt)},
qual_ {qual},
filter_ {std::forward<Filters>(filters)},
info_ {std::forward<Info>(info)},
format_ {},
genotypes_ {},
samples_ {}
{}

template <typename String1, typename String2, typename Sequence1, typename Sequence2,
          typename Filters, typename Info, typename Format, typename Genotypes, typename Samples>
VcfRecord::VcfRecord(String1&& chrom, SizeType pos, String2&& id,   Sequence1&& ref, Sequence2&& alt,
                     boost::optional<QualityType> qual, Filters&& filters,
                     Info&& info, Format&& format, Genotypes&& genotypes, Samples&& samples)
:
chrom_ {std::forward<String1>(chrom)},
pos_ {pos},
id_ {std::forward<String2>(id)},
ref_ {std::forward<Sequence1>(ref)},
alt_ {std::forward<Sequence2>(alt)},
qual_ {qual},
filter_ {std::forward<Filters>(filters)},
info_ {std::forward<Info>(info)},
format_ {std::forward<Format>(format)},
genotypes_ {std::forward<Genotypes>(genotypes)},
samples_ {std::forward<Samples>(samples)}
{}

// non-member functions

VcfRecord::NucleotideSequence get_ancestral_allele(const VcfRecord& record);
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

std::ostream& operator<<(std::ostream& os, const VcfRecord& record);

class VcfRecord::Builder
{
public:
    using SizeType           = VcfRecord::SizeType;
    using IdType             = VcfRecord::IdType;
    using NucleotideSequence = VcfRecord::NucleotideSequence;
    using QualityType        = VcfRecord::QualityType;
    using SampleName       = VcfRecord::SampleName;
    using KeyType            = VcfRecord::KeyType;
    using ValueType          = VcfRecord::ValueType;
    
    enum class Phasing { Phased, Unphased };
    
    Builder() = default;
    
    Builder(const VcfRecord& call);
    
    Builder& set_chrom(std::string name);
    Builder& set_pos(SizeType pos);
    Builder& set_id(IdType id);
    Builder& set_ref(const char allele);
    Builder& set_ref(NucleotideSequence allele);
    Builder& set_alt(const char allele); // if just one
    Builder& set_alt(NucleotideSequence allele); // if just one
    Builder& set_alt(std::vector<NucleotideSequence> alleles);
    Builder& set_qual(QualityType quality);
    Builder& set_passed();
    Builder& set_filter(std::vector<KeyType> filter);
    Builder& set_filter(std::initializer_list<KeyType> filter);
    Builder& add_filter(KeyType filter);
    Builder& reserve_info(unsigned n);
    Builder& add_info(const KeyType& key); // flags
    Builder& set_info(const KeyType& key, const ValueType& value);
    template <typename T> Builder& set_info(const KeyType& key, const T& value); // calls std::to_string
    Builder& set_info(const KeyType& key, std::vector<ValueType> values);
    Builder& set_info(const KeyType& key, std::initializer_list<ValueType> values);
    Builder& set_info_flag(KeyType key);
    Builder& clear_info() noexcept;
    Builder& clear_info(const KeyType& key);
    Builder& set_format(std::vector<KeyType> format);
    Builder& set_format(std::initializer_list<KeyType> format);
    Builder& add_format(KeyType key);
    Builder& set_homozygous_ref_genotype(const SampleName& sample, unsigned ploidy);
    Builder& reserve_samples(unsigned n);
    Builder& set_genotype(const SampleName& sample, std::vector<NucleotideSequence> alleles, Phasing phasing);
    Builder& set_genotype(const SampleName& sample, const std::vector<boost::optional<unsigned>>& alleles, Phasing is_phased);
    Builder& set_format(const SampleName& sample, const KeyType& key, const ValueType& value);
    template <typename T>
    Builder& set_format(const SampleName& sample, const KeyType& key, const T& value); // calls std::to_string
    Builder& set_format(const SampleName& sample, const KeyType& key, std::vector<ValueType> values);
    Builder& set_format(const SampleName& sample, const KeyType& key, std::initializer_list<ValueType> values);
    Builder& set_format_missing(const SampleName& sample, const KeyType& key);
    
    Builder& set_refcall();
    Builder& set_somatic();
    
    SizeType pos() const noexcept;
    
    VcfRecord build() const;
    VcfRecord build_once() noexcept;
    
private:
    decltype(VcfRecord::chrom_) chrom_ = ".";
    decltype(VcfRecord::pos_) pos_ = 0;
    decltype(VcfRecord::id_) id_ = ".";
    decltype(VcfRecord::ref_) ref_ = ".";
    decltype(VcfRecord::alt_) alt_ = {"."};
    decltype(VcfRecord::qual_) qual_ = boost::none;
    decltype(VcfRecord::filter_) filter_ = {};
    decltype(VcfRecord::info_) info_ = {};
    decltype(VcfRecord::format_) format_ = {};
    decltype(VcfRecord::genotypes_) genotypes_ = {};
    decltype(VcfRecord::samples_) samples_ = {};
};

template <typename T>
VcfRecord::Builder& VcfRecord::Builder::set_info(const KeyType& key, const T& value)
{
    return set_info(key, std::to_string(value));
}

template <typename T>
VcfRecord::Builder& VcfRecord::Builder::set_format(const SampleName& sample, const KeyType& key,
                                                   const T& value)
{
    return set_format(sample, key, std::to_string(value));
}

} // namespace octopus

namespace std {
    template <> struct hash<octopus::VcfRecord>
    {
        size_t operator()(const octopus::VcfRecord& record) const
        {
            return hash<string>()(record.id());
        }
    };
} // namespace std

#endif /* defined(__Octopus__vcf_record__) */
