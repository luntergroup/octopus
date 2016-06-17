//
//  read_transformations.hpp
//  Octopus
//
//  Created by Daniel Cooke on 07/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_read_transformations_hpp
#define Octopus_read_transformations_hpp

#include "aligned_read.hpp"

namespace Octopus { namespace ReadTransforms
{
    struct CapBaseQualities
    {
        using QualityType = AlignedRead::QualityType;
        
        CapBaseQualities() = default;
        
        explicit CapBaseQualities(QualityType max_quality);
        
        void operator()(AlignedRead& read) const noexcept;
        
    private:
        QualityType max_quality_;
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
        using SizeType = AlignedRead::SizeType;
        
        MaskTail() = default;
        
        explicit MaskTail(SizeType num_bases);
        
        void operator()(AlignedRead& read) const noexcept;
        
    private:
        const SizeType num_bases_;
    };
    
    struct MaskSoftClipped
    {
        void operator()(AlignedRead& read) const noexcept;
    };
    
    struct MaskSoftClippedBoundries
    {
        using SizeType = AlignedRead::SizeType;
        
        MaskSoftClippedBoundries() = default;
        
        explicit MaskSoftClippedBoundries(SizeType num_bases);
        
        void operator()(AlignedRead& read) const noexcept;
        
    private:
        const SizeType num_bases_;
    };
} // namespace ReadTransforms
} // namespace Octopus

#endif
