// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_transform_hpp
#define read_transform_hpp

#include <functional>
#include <vector>

#include "basics/aligned_read.hpp"
#include "io/reference/reference_genome.hpp"

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
    Length num_bases_;
};

struct MaskLowQualityTails
{
    using BaseQuality = AlignedRead::BaseQuality;
    
    MaskLowQualityTails() = default;
    
    explicit MaskLowQualityTails(BaseQuality threshold);
    
    void operator()(AlignedRead& read) const noexcept;

private:
    BaseQuality threshold_;
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
    Length num_bases_;
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

struct MaskLowQualitySoftClippedBoundaryBases
{
    using Length = AlignedRead::NucleotideSequence::size_type;
    using BaseQuality = AlignedRead::BaseQuality;
    
    MaskLowQualitySoftClippedBoundaryBases() = default;
    
    explicit MaskLowQualitySoftClippedBoundaryBases(Length num_bases, BaseQuality max);
    
    void operator()(AlignedRead& read) const noexcept;

private:
    Length num_bases_;
    BaseQuality max_;
};

struct MaskLowAverageQualitySoftClippedTails
{
    using BaseQuality = AlignedRead::BaseQuality;
    using Length = AlignedRead::NucleotideSequence::size_type;
    
    MaskLowAverageQualitySoftClippedTails() = default;
    
    explicit MaskLowAverageQualitySoftClippedTails(BaseQuality threshold, Length min_tail_length = 1);
    
    void operator()(AlignedRead& read) const noexcept;

private:
    BaseQuality threshold_;
    Length min_tail_length_;
};

struct MaskInvertedSoftClippedReadEnds
{
    MaskInvertedSoftClippedReadEnds(const ReferenceGenome& reference,
                                    AlignedRead::NucleotideSequence::size_type min_clip_length,
                                    GenomicRegion::Size max_flank_search);
    
    void operator()(AlignedRead& read) const;
    
private:
    std::reference_wrapper<const ReferenceGenome> reference_;
    AlignedRead::NucleotideSequence::size_type min_clip_length_;
    GenomicRegion::Size max_flank_search_;
};

struct Mask3PrimeShiftedSoftClippedHeads
{
    Mask3PrimeShiftedSoftClippedHeads(const ReferenceGenome& reference,
                                      AlignedRead::NucleotideSequence::size_type min_clip_length,
                                      GenomicRegion::Size max_flank_search);
    
    void operator()(AlignedRead& read) const;

private:
    std::reference_wrapper<const ReferenceGenome> reference_;
    AlignedRead::NucleotideSequence::size_type min_clip_length_;
    GenomicRegion::Size max_flank_search_;
};

struct ClearAnnotations
{
    void operator()(AlignedRead& read) const noexcept;
};

using ReadReferenceVector = std::vector<std::reference_wrapper<AlignedRead>>;

struct MaskTemplateAdapters
{
    void operator()(ReadReferenceVector& read_template) const;
};

struct MaskStrandOfDuplicatedBases
{
    void operator()(ReadReferenceVector& read_template) const;
};

struct MaskClippedDuplicatedBases
{
    void operator()(ReadReferenceVector& read_template) const;
};

} // namespace readpipe
} // namespace octopus

#endif
