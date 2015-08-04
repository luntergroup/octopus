//
//  vcf_record.h
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
#include <ostream>

class VcfRecord
{
public:
    using SizeType     = std::uint_fast32_t;
    using SequenceType = std::string;
    using QualityType  = std::uint_fast8_t;
    
    VcfRecord()  = default;
    
    // constructor without genotype fields
    template <typename StringType1_, typename StringType2_, typename SequenceType1_, typename SequenceType2_,
              typename Filters_, typename Info_>
    VcfRecord(StringType1_&& chrom, SizeType pos, StringType2_&& id, SequenceType1_&& ref, SequenceType2_&& alt,
              QualityType qual, Filters_&& filters, Info_&& info);
    
    // constructor with genotype fields
    template <typename StringType1_, typename StringType2_, typename SequenceType1_, typename SequenceType2_,
    typename Filters_, typename Info_, typename Format_, typename Genotypes_>
    VcfRecord(StringType1_&& chrom, SizeType pos, StringType2_&& id, SequenceType1_&& ref, SequenceType2_&& alt,
              QualityType qual, Filters_&& filters, Info_&& info, Format_&& format, Genotypes_&& genotypes);
    
    ~VcfRecord() = default;
    
    VcfRecord(const VcfRecord&)            = default;
    VcfRecord& operator=(const VcfRecord&) = default;
    VcfRecord(VcfRecord&&)                 = default;
    VcfRecord& operator=(VcfRecord&&)      = default;
    
    const std::string& get_chromosome_name() const noexcept;
    SizeType get_position() const noexcept;
    const std::string& get_id() const noexcept;
    const SequenceType& get_ref_allele() const noexcept;
    unsigned get_num_alt_alleles() const noexcept;
    const SequenceType& get_alt_allele(unsigned n) const noexcept;
    
    QualityType get_quality() const noexcept;
    bool has_filter(const std::string& filter) const noexcept;
    unsigned get_num_samples() const noexcept;
    
    bool has_info(const std::string& key) const noexcept;
    const std::vector<std::string>& get_info_values(const std::string& key) const;
    
    bool has_genotype_data() const noexcept;
    
    friend std::ostream& operator<<(std::ostream& os, const VcfRecord& record);
    
private:
    // mandatory fields
    std::string chrom_;
    SizeType pos_;
    std::string id_;
    SequenceType ref_;
    std::vector<SequenceType> alt_;
    QualityType qual_;
    std::vector<std::string> filter_;
    std::map<std::string, std::vector<std::string>> info_;
    
    // optional fields
    std::vector<std::string> format_;
    std::map<std::string, std::vector<std::string>> genotypes_;
};

template <typename StringType1_, typename StringType2_, typename SequenceType1_, typename SequenceType2_,
          typename Filters_, typename Info_>
VcfRecord::VcfRecord(StringType1_&& chrom, SizeType pos, StringType2_&& id, SequenceType1_&& ref,
                     SequenceType2_&& alt, QualityType qual, Filters_&& filters, Info_&& info)
:
chrom_ {std::forward<StringType1_>(chrom)},
pos_ {pos},
id_ {std::forward<StringType2_>(id)},
ref_ {std::forward<SequenceType1_>(ref)},
alt_ {std::forward<SequenceType2_>(alt)},
qual_ {qual},
filter_ {std::forward<Filters_>(filters)},
info_ {std::forward<Info_>(info)},
format_ {},
genotypes_ {}
{}

template <typename StringType1_, typename StringType2_, typename SequenceType1_, typename SequenceType2_,
          typename Filters_, typename Info_, typename Format_, typename Genotypes_>
VcfRecord::VcfRecord(StringType1_&& chrom, SizeType pos, StringType2_&& id, SequenceType1_&& ref, SequenceType2_&& alt,
                     QualityType qual, Filters_&& filters, Info_&& info, Format_&& format, Genotypes_&& genotypes)
:
chrom_ {std::forward<StringType1_>(chrom)},
pos_ {pos},
id_ {std::forward<StringType2_>(id)},
ref_ {std::forward<SequenceType1_>(ref)},
alt_ {std::forward<SequenceType2_>(alt)},
qual_ {qual},
filter_ {std::forward<Filters_>(filters)},
info_ {std::forward<Info_>(info)},
format_ {std::forward<Format_>(format)},
genotypes_ {std::forward<Genotypes_>(genotypes)}
{}

// non-member functions

bool is_dbsnp_member(const VcfRecord& record) noexcept;
bool is_hapmap2_member(const VcfRecord& record) noexcept;
bool is_hapmap3_member(const VcfRecord& record) noexcept;
bool is_1000g_member(const VcfRecord& record) noexcept;
bool is_somatic(const VcfRecord& record) noexcept;
bool is_validated(const VcfRecord& record) noexcept;

std::ostream& operator<<(std::ostream& os, const VcfRecord& record);

#endif /* defined(__Octopus__vcf_record__) */
