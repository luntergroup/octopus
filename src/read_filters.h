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
#include <iterator>  // std::prev, std::cbegin, std::cend

#include "aligned_read.h"
#include "cigar_string.h"

// Context-free filters

inline bool is_not_secondary_alignment(const AlignedRead& read)
{
    return !read.is_marked_secondary_alignment();
}

inline bool is_good_mapping_quality(const AlignedRead& read,
                                    AlignedRead::QualityType min_mapping_quality)
{
    return read.get_mapping_quality() >= min_mapping_quality;
}

inline bool has_sufficient_good_quality_bases(const AlignedRead& read,
                                              AlignedRead::QualityType min_base_quality,
                                              unsigned min_num_good_quality_bases)
{
    const auto& qualities = read.get_qualities();
    return std::count_if(std::cbegin(qualities), std::cend(qualities), [min_base_quality] (const auto& qual) {
        return qual >= min_base_quality;
    }) >= min_num_good_quality_bases;
}

inline bool is_mapped(const AlignedRead& read)
{
    return !read.is_marked_unmapped();
}

inline bool is_not_chimeric(const AlignedRead& read)
{
    return !read.is_chimeric();
}

//inline bool is_next_segment_mapped(const AlignedRead& read)
//{
//    return (read.has_mate_pair()) ? read.is_marked_proper_pair() : true;
//}

inline bool is_not_marked_duplicate(const AlignedRead& read)
{
    return !read.is_marked_duplicate();
}

inline bool is_short(const AlignedRead& read, AlignedRead::SizeType max_length)
{
    return read.get_sequence_size() <= max_length;
}

inline bool is_long(const AlignedRead& read, AlignedRead::SizeType min_length)
{
    return read.get_sequence_size() <= min_length;
}

inline bool is_not_contaminated(const AlignedRead& read)
{
    return (read.is_chimeric()) ?
        read.get_sequence_size() >= read.get_next_segment()->get_inferred_template_length() : true;
}

inline bool not_marked_qc_fail(const AlignedRead& read)
{
    return !read.is_marked_qc_fail();
}

// Context-based filters

template <typename BidirectionalIterator>
bool is_not_duplicate(const AlignedRead& read, BidirectionalIterator first_good_read,
                      BidirectionalIterator previous_good_read)
{
    if (first_good_read != previous_good_read) {
        if (read.is_chimeric() && previous_good_read->is_chimeric()) {
            return !(read == *previous_good_read &&
                    *read.get_next_segment() == *(previous_good_read->get_next_segment()));
        } else {
            return read != *previous_good_read;
        }
    } else {
        return true;
    }
}

#endif /* defined(__Octopus__read_filters__) */
