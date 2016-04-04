//
//  read_transformations.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "read_transformations.hpp"

#include <iostream>

namespace Octopus
{
namespace ReadTransforms
{
    void trim_overlapping::operator()(AlignedRead& read) const
    {
        if (read.is_chimeric() && !read.is_marked_reverse_mapped()
            && read.get_next_segment().get_begin() < region_end(read)) {
            const auto overlapped_size = region_end(read) - read.get_next_segment().get_begin();
            read.zero_back_qualities(overlapped_size);
        }
    }
    
    void trim_adapters::operator()(AlignedRead& read) const
    {
        if (read.is_chimeric()) {
            const auto insert_size = read.get_next_segment().get_inferred_template_length();
            const auto read_size   = sequence_size(read);
            
            if (insert_size <= read_size) {
                const auto num_adapter_bases = read_size - insert_size;
                
                if (read.is_marked_reverse_mapped()) {
                    read.zero_back_qualities(num_adapter_bases);
                } else {
                    read.zero_front_qualities(num_adapter_bases);
                }
            }
        }
    }
    
    trim_tail::trim_tail(SizeType num_bases) : num_bases_ {num_bases} {};
    
    void trim_tail::operator()(AlignedRead& read) const
    {
        if (read.is_marked_reverse_mapped()) {
            read.zero_front_qualities(num_bases_);
        } else {
            read.zero_back_qualities(num_bases_);
        }
    }
    
    void trim_soft_clipped::operator()(AlignedRead& read) const noexcept
    {
        if (is_soft_clipped(read.get_cigar_string())) {
            const auto soft_clipped_sizes = get_soft_clipped_sizes(read.get_cigar_string());
            read.zero_front_qualities(soft_clipped_sizes.first);
            read.zero_back_qualities(soft_clipped_sizes.second);
        }
    }
    
    trim_soft_clipped_tails::trim_soft_clipped_tails(SizeType num_bases) : num_bases_ {num_bases} {};
    
    void trim_soft_clipped_tails::operator()(AlignedRead& read) const
    {
        if (is_soft_clipped(read)) {
            const auto soft_clipped_sizes = get_soft_clipped_sizes(read);
            read.zero_front_qualities(soft_clipped_sizes.first + num_bases_);
            read.zero_back_qualities(soft_clipped_sizes.second + num_bases_);
        } else {
            trim_tail{num_bases_}(read);
        }
    }
} // namespace ReadTransforms
} // namespace Octopus
