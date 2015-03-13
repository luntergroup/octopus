//
//  read_filters.h
//  Octopus
//
//  Created by Daniel Cooke on 07/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_filters__
#define __Octopus__read_filters__

#include <functional>
#include <algorithm> // std::count_if
#include <iterator>  // std::prev

#include "aligned_read.h"
#include "cigar_string.h"

using std::cbegin;
using std::cend;

// Context-free filters

inline bool is_not_secondary_alignment(const AlignedRead& a_read)
{
    return !a_read.is_marked_secondary_alignment();
}

inline bool is_good_mapping_quality(const AlignedRead& a_read,
                                    AlignedRead::QualityType min_mapping_quality)
{
    return a_read.get_mapping_quality() >= min_mapping_quality;
}

inline bool has_sufficient_good_quality_bases(const AlignedRead& a_read,
                                              AlignedRead::QualityType min_base_quality,
                                              unsigned min_num_good_quality_bases)
{
    const auto& qualities = a_read.get_qualities();
    return std::count_if(cbegin(qualities), cend(qualities), [min_base_quality] (const auto& qual) {
        return qual >= min_base_quality;
    }) >= min_num_good_quality_bases;
}

inline bool is_mapped(const AlignedRead& a_read)
{
    return !a_read.is_marked_unmapped();
}

inline bool is_not_chimeric(const AlignedRead& a_read)
{
    return !a_read.is_chimeric();
}

//inline bool is_next_segment_mapped(const AlignedRead& a_read)
//{
//    return (a_read.has_mate_pair()) ? a_read.is_marked_proper_pair() : true;
//}

inline bool is_not_marked_duplicate(const AlignedRead& a_read)
{
    return !a_read.is_marked_duplicate();
}

inline bool is_not_contaminated(const AlignedRead& a_read)
{
    return (a_read.is_chimeric()) ?
        a_read.get_sequence_size() >= a_read.get_next_segment()->get_inferred_template_length() : true;
}

// Context-based filters

template <typename BidirectionalIterator>
bool is_not_duplicate(const AlignedRead& a_read, BidirectionalIterator first_good_read,
                      BidirectionalIterator previous_good_read)
{
    if (first_good_read != previous_good_read) {
        if (a_read.is_chimeric() && previous_good_read->is_chimeric()) {
            return !(a_read == *previous_good_read &&
                    *a_read.get_next_segment() == *(previous_good_read->get_next_segment()));
        } else {
            return a_read != *previous_good_read;
        }
    } else {
        return true;
    }
}

#endif /* defined(__Octopus__read_filters__) */
