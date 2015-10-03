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

namespace Octopus { namespace ReadTransforms {
    
    struct trim_adapters
    {
        trim_adapters() = default;
        
        void operator()(AlignedRead& read) const
        {
            if (read.is_chimeric()) {
                auto insert_size = read.get_next_segment().get_inferred_template_length();
                auto read_size   = read.get_sequence_size();
                
                if (insert_size > 0 && insert_size < read_size) {
                    auto num_adapter_bases = read_size - insert_size;
                    if (read.is_marked_reverse_mapped()) {
                        read.zero_back_qualities(num_adapter_bases);
                    } else {
                        read.zero_front_qualities(num_adapter_bases);
                    }
                }
            }
        }
    };
    
    // TODO: perhaps it would be better to just set these bases to reference rather
    // than zeroing... this could avoid a lot of spurious candidates. Would need to add
    // sequence setting api to AlignedRead
    struct trim_tail
    {
        using SizeType = AlignedRead::SizeType;
        
        trim_tail() = default;
        explicit trim_tail(SizeType num_bases) : num_bases_ {num_bases} {};
        
        void operator()(AlignedRead& read) const
        {
            if (read.is_marked_reverse_mapped()) {
                read.zero_front_qualities(num_bases_);
            } else {
                read.zero_back_qualities(num_bases_);
            }
        }
        
    private:
        const SizeType num_bases_;
    };
    
    struct trim_soft_clipped
    {
        trim_soft_clipped() = default;
        
        void operator()(AlignedRead& read) const
        {
            if (is_soft_clipped(read.get_cigar_string())) {
                auto qualities = read.get_qualities();
                auto soft_clipped_sizes = get_soft_clipped_sizes(read.get_cigar_string());
                read.zero_front_qualities(soft_clipped_sizes.first);
                read.zero_back_qualities(soft_clipped_sizes.second);
            }
        }
    };

} // namespace ReadTransforms
} // namespace Octopus

#endif
