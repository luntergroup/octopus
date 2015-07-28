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

enum class Filter { PASS, q10};

class VcfRecord
{
public:
    using SizeType     = std::uint_fast32_t;
    using SequenceType = std::string;
    using QualityType  = std::uint_fast8_t;
    
    VcfRecord()  = default;
    ~VcfRecord() = default;
    
    VcfRecord(const VcfRecord&)            = default;
    VcfRecord& operator=(const VcfRecord&) = default;
    VcfRecord(VcfRecord&&)                 = default;
    VcfRecord& operator=(VcfRecord&&)      = default;
    
private:
    // mandatory fields
    std::string chrom_;
    SizeType pos_;
    std::string id_;
    SequenceType ref_;
    SequenceType alt_;
    QualityType qual_;
    std::string filter; // should be a string?
    //std::vector<>
    
    // optional fields
};

#endif /* defined(__Octopus__vcf_record__) */
