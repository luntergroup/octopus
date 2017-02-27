// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_transform_hpp
#define read_transform_hpp

#include <functional>
#include <vector>

#include "basics/aligned_read.hpp"

namespace octopus { namespace readpipe {

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
    const BaseQuality max_;
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

struct MaskLowQualityTails
{
    using BaseQuality = AlignedRead::BaseQuality;
    
    MaskLowQualityTails() = default;
    
    explicit MaskLowQualityTails(BaseQuality threshold);
    
    void operator()(AlignedRead& read) const noexcept;

private:
    const BaseQuality threshold_;
};

struct MaskSoftClipped
{
    void operator()(AlignedRead& read) const noexcept;
};

struct MaskSoftClippedBoundraryBases
{
    using Length = AlignedRead::NucleotideSequence::size_type;
    
    MaskSoftClippedBoundraryBases() = default;
    
    explicit MaskSoftClippedBoundraryBases(Length num_bases);
    
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
    const BaseQuality max_;
};

struct MaskLowQualitySoftClippedBoundaryBases
{
    using Length = AlignedRead::NucleotideSequence::size_type;
    using BaseQuality = AlignedRead::BaseQuality;
    
    MaskLowQualitySoftClippedBoundaryBases() = default;
    
    explicit MaskLowQualitySoftClippedBoundaryBases(Length num_bases, BaseQuality max);
    
    void operator()(AlignedRead& read) const noexcept;

private:
    const Length num_bases_;
    const BaseQuality max_;
};

struct MaskLowAverageQualitySoftClippedTails
{
    using BaseQuality = AlignedRead::BaseQuality;
    using Length = AlignedRead::NucleotideSequence::size_type;
    
    MaskLowAverageQualitySoftClippedTails() = default;
    
    explicit MaskLowAverageQualitySoftClippedTails(BaseQuality threshold, Length min_tail_length = 1);
    
    void operator()(AlignedRead& read) const noexcept;

private:
    const BaseQuality threshold_;
    const Length min_tail_length_;
};

using ReadReferenceVector = std::vector<std::reference_wrapper<AlignedRead>>;

struct MaskTemplateAdapters
{
    void operator()(ReadReferenceVector& read_template) const;
};

struct MaskDuplicatedBases
{
    void operator()(ReadReferenceVector& read_template) const;
};

} // namespace readpipe
} // namespace octopus

#endif
