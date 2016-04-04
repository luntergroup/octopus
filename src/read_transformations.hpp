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
    struct trim_overlapping
    {
        void operator()(AlignedRead& read) const;
    };
    
    struct trim_adapters
    {
        void operator()(AlignedRead& read) const;
    };
    
    struct trim_tail
    {
        using SizeType = AlignedRead::SizeType;
        
        trim_tail() = default;
        explicit trim_tail(SizeType num_bases);
        
        void operator()(AlignedRead& read) const;
        
    private:
        const SizeType num_bases_;
    };
    
    struct trim_soft_clipped
    {
        void operator()(AlignedRead& read) const noexcept;
    };
    
    struct trim_soft_clipped_tails
    {
        using SizeType = AlignedRead::SizeType;
        
        trim_soft_clipped_tails() = default;
        explicit trim_soft_clipped_tails(SizeType num_bases);
        
        void operator()(AlignedRead& read) const;
        
    private:
        const SizeType num_bases_;
    };
} // namespace ReadTransforms
} // namespace Octopus

#endif
