//
//  aligned_read_builder.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "aligned_read_builder.hpp"

#include <utility>

AlignedReadBuilder& AlignedReadBuilder::set_region(GenomicRegion region)
{
    region_ = std::move(region);
    return *this;
}

AlignedReadBuilder& AlignedReadBuilder::set_sequence(SequenceType sequence)
{
    sequence_ = std::move(sequence);
    return *this;
}

AlignedReadBuilder& AlignedReadBuilder::set_qualities(Qualities qualities)
{
    qualities_ = std::move(qualities);
    return *this;
}

AlignedReadBuilder& AlignedReadBuilder::set_cigar(CigarString cigar)
{
    cigar_string_ = std::move(cigar);
    return *this;
}

AlignedReadBuilder& AlignedReadBuilder::set_read_group(ReadGroupType read_group)
{
    read_group_ = std::move(read_group);
    return *this;
}

AlignedReadBuilder& AlignedReadBuilder::set_flags(Flags flags)
{
    flags_ = std::move(flags);
    return *this;
}

AlignedReadBuilder& AlignedReadBuilder::set_mapping_quality(QualityType quality)
{
    mapping_quality_ = std::move(quality);
    return *this;
}

AlignedRead AlignedReadBuilder::build() const
{
    return AlignedRead
    {
        region_,
        sequence_,
        qualities_,
        cigar_string_,
        mapping_quality_,
        flags_
    };
}

AlignedRead AlignedReadBuilder::build_once()
{
    return AlignedRead
    {
        std::move(region_),
        std::move(sequence_),
        std::move(qualities_),
        std::move(cigar_string_),
        std::move(mapping_quality_),
        std::move(flags_)
    };
}
