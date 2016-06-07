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
    struct MaskOverlappedSegment
    {
        void operator()(AlignedRead& read) const noexcept;
    };
    
    struct MaskAdapters
    {
        void operator()(AlignedRead& read) const noexcept;
    };
    
    struct TrimTail
    {
        using SizeType = AlignedRead::SizeType;
        
        TrimTail() = default;
        explicit TrimTail(SizeType num_bases);
        
        void operator()(AlignedRead& read) const noexcept;
        
    private:
        const SizeType num_bases_;
    };
    
    struct MaskSoftClipped
    {
        void operator()(AlignedRead& read) const noexcept;
    };
    
    struct TrimSoftClippedTails
    {
        using SizeType = AlignedRead::SizeType;
        
        TrimSoftClippedTails() = default;
        explicit TrimSoftClippedTails(SizeType num_bases);
        
        void operator()(AlignedRead& read) const noexcept;
        
    private:
        const SizeType num_bases_;
    };
} // namespace ReadTransforms
} // namespace Octopus

#endif
