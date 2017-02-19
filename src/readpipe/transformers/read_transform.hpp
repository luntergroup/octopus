// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_transform_hpp
#define read_transform_hpp

#include "basics/aligned_read.hpp"

namespace octopus { namespace readpipe
{

struct CapitaliseBases
{
    void operator()(AlignedRead& read) const noexcept;
};

struct CapBaseQualities
{
    using BaseQuality = AlignedRead::BaseQuality;
    
    CapBaseQualities() = default;
    
    explicit CapBaseQualities(BaseQuality max);
    
    void operator()(AlignedRead& read) const noexcept;
    
private:
    BaseQuality max_;
};

struct MaskOverlappedSegment
{
    void operator()(AlignedRead& read) const noexcept;
};

struct MaskAdapters
{
    void operator()(AlignedRead& read) const noexcept;
};

struct MaskTail
{
    using Length = AlignedRead::NucleotideSequence::size_type;
    
    MaskTail() = default;
    
    explicit MaskTail(Length num_bases);
    
    void operator()(AlignedRead& read) const noexcept;
    
private:
    const Length num_bases_;
};

struct MaskSoftClipped
{
    void operator()(AlignedRead& read) const noexcept;
};

struct MaskSoftClippedBoundries
{
    using Length = AlignedRead::NucleotideSequence::size_type;
    
    MaskSoftClippedBoundries() = default;
    
    explicit MaskSoftClippedBoundries(Length num_bases);
    
    void operator()(AlignedRead& read) const noexcept;
    
private:
    const Length num_bases_;
};

struct MaskLowQualitySoftClippedBases
{
    using BaseQuality = AlignedRead::BaseQuality;
    
    MaskLowQualitySoftClippedBases() = default;
    
    explicit MaskLowQualitySoftClippedBases(BaseQuality max);
    
    void operator()(AlignedRead& read) const noexcept;

private:
    BaseQuality max_;
};

struct QualityAdjustedSoftClippedMasker
{
    void operator()(AlignedRead& read) const noexcept;
};

} // namespace readpipe
} // namespace octopus

#endif
