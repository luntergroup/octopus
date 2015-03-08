//
//  read_filters.h
//  Octopus
//
//  Created by Daniel Cooke on 07/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include <functional>
#include <algorithm> // std::count_if
#include <iterator>  // std::prev

#include "aligned_read.h"
#include "cigar_string.h"

using std::cbegin;
using std::cend;

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

inline bool is_good_pair(const AlignedRead& a_read)
{
    return (a_read.has_mate_pair()) ? a_read.is_marked_proper_pair() : true;
}

inline bool is_not_marked_duplicate(const AlignedRead& a_read)
{
    return !a_read.is_marked_duplicate();
}

template <typename BidirectionalIterator>
bool is_not_duplicate(const AlignedRead& a_read, BidirectionalIterator first_good,
                      BidirectionalIterator previous_good)
{
    if (first_good != previous_good) {
        if (a_read.has_mate_pair() && previous_good->has_mate_pair()) {
            return !(a_read == *previous_good &&
                    *a_read.get_mate_pair() == *(previous_good->get_mate_pair()));
        } else {
            return a_read != *previous_good;
        }
    } else {
        return true;
    }
}
