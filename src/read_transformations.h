//
//  read_transformations.h
//  Octopus
//
//  Created by Daniel Cooke on 07/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_read_transformations_h
#define Octopus_read_transformations_h

#include "aligned_read.h"

inline void trim_adapters(AlignedRead& a_read)
{
    if (a_read.is_chimeric()) {
        auto insert_size = a_read.get_next_segment()->get_inferred_template_length();
        auto read_size   = a_read.get_sequence_size();
        if (insert_size > 0 && insert_size < read_size) {
            auto num_adapter_bases = read_size - insert_size;
            if (a_read.is_marked_reverse_mapped()) {
                a_read.zero_back_qualities(num_adapter_bases);
            } else {
                a_read.zero_front_qualities(num_adapter_bases);
            }
        }
    }
}

inline void trim_soft_clipped(AlignedRead& a_read)
{
    if (is_soft_clipped(a_read.get_cigar_string())) {
        auto qualities = a_read.get_qualities();
        auto soft_clipped_sizes = get_soft_clipped_sizes(a_read.get_cigar_string());
        a_read.zero_front_qualities(soft_clipped_sizes.first);
        a_read.zero_back_qualities(soft_clipped_sizes.second);
    }
}

#endif
