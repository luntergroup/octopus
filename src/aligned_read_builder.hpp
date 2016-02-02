//
//  aligned_read_builder.hpp
//  Octopus
//
//  Created by Daniel Cooke on 01/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef aligned_read_builder_hpp
#define aligned_read_builder_hpp

#include "aligned_read.hpp"
#include "genomic_region.hpp"
#include "cigar_string.hpp"

class AlignedReadBuilder
{
public:
    using SizeType      = AlignedRead::SizeType;
    using SequenceType  = AlignedRead::SequenceType;
    using QualityType   = AlignedRead::QualityType;
    using Qualities     = AlignedRead::Qualities;
    using ReadGroupType = AlignedRead::ReadGroupType;
    using Flags         = AlignedRead::Flags;
    using NextSegment   = AlignedRead::NextSegment;
    
    AlignedReadBuilder()  = default;
    ~AlignedReadBuilder() = default;
    
    AlignedReadBuilder(const AlignedReadBuilder&)            = default;
    AlignedReadBuilder& operator=(const AlignedReadBuilder&) = default;
    AlignedReadBuilder(AlignedReadBuilder&&)                 = default;
    AlignedReadBuilder& operator=(AlignedReadBuilder&&)      = default;
    
    AlignedReadBuilder& set_region(GenomicRegion region);
    AlignedReadBuilder& set_sequence(SequenceType sequence);
    AlignedReadBuilder& set_qualities(Qualities qualities);
    AlignedReadBuilder& set_cigar(CigarString cigar);
    AlignedReadBuilder& set_read_group(ReadGroupType read_group);
    AlignedReadBuilder& set_flags(Flags flags);
    AlignedReadBuilder& set_mapping_quality(QualityType quality);
    
    AlignedRead build() const;
    AlignedRead build_once();
private:
    GenomicRegion region_;
    SequenceType sequence_ = "";
    Qualities qualities_ = {};
    CigarString cigar_string_ = {};
    ReadGroupType read_group_ = "";
    //std::unique_ptr<NextSegment> next_segment_;
    Flags flags_;
    QualityType mapping_quality_ = 0;
};

#endif /* aligned_read_builder_hpp */
