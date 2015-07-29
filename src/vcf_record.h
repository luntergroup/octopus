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
#include <ostream>

class VcfRecord
{
public:
    using SizeType     = std::uint_fast32_t;
    using SequenceType = std::string;
    using QualityType  = std::uint_fast8_t;
    
    VcfRecord()  = default;
    template <typename StringType1_, typename StringType2_, typename SequenceType1_, typename SequenceType2_>
    VcfRecord(StringType1_&& chrom, SizeType pos, StringType2_&& id, SequenceType1_&& ref, SequenceType2_&& alt, QualityType qual);
    ~VcfRecord() = default;
    
    VcfRecord(const VcfRecord&)            = default;
    VcfRecord& operator=(const VcfRecord&) = default;
    VcfRecord(VcfRecord&&)                 = default;
    VcfRecord& operator=(VcfRecord&&)      = default;
    
    const std::string& get_chrom() const noexcept;
    SizeType get_pos() const noexcept;
    const std::string& get_id() const noexcept;
    const SequenceType& get_ref() const noexcept;
    unsigned get_num_alt_alleles() const noexcept;
    const SequenceType& get_alt(unsigned n) const noexcept;
    
    QualityType get_qual() const noexcept;
    bool has_filter(const std::string& filter) const noexcept;
    unsigned get_num_samples() const noexcept;
    
private:
    // mandatory fields
    std::string chrom_;
    SizeType pos_;
    std::string id_;
    SequenceType ref_;
    std::vector<SequenceType> alt_;
    QualityType qual_;
    std::vector<std::string> filter_;
    std::vector<std::string> info_;
    
    // optional fields
    std::vector<std::string> format_;
    std::vector<std::vector<std::string>> genotypes_;
};

template <typename StringType1_, typename StringType2_, typename SequenceType1_, typename SequenceType2_>
VcfRecord::VcfRecord(StringType1_&& chrom, SizeType pos, StringType2_&& id, SequenceType1_&& ref,
                     SequenceType2_&& alt, QualityType qual)
:
chrom_ {std::forward<StringType1_>(chrom)},
pos_ {pos},
id_ {std::forward<StringType2_>(id)},
ref_ {std::forward<SequenceType1_>(ref)},
alt_ {std::forward<SequenceType2_>(alt)},
qual_ {qual}
{}

std::ostream& operator<<(std::ostream& os, const VcfRecord& record);

#endif /* defined(__Octopus__vcf_record__) */
