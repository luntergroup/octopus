//
//  read_filters.h
//  Octopus
//
//  Created by Daniel Cooke on 07/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_filters__
#define __Octopus__read_filters__

#include <algorithm>
#include <iterator>

#include "aligned_read.hpp"
#include "cigar_string.hpp"

namespace Octopus { namespace ReadFilters {

// Context-free filters

struct is_not_secondary_alignment
{
    is_not_secondary_alignment() = default;
    
    bool operator()(const AlignedRead& read) const
    {
        return !read.is_marked_secondary_alignment();
    }
};
    
struct is_not_supplementary_alignment
{
    is_not_supplementary_alignment() = default;
    
    bool operator()(const AlignedRead& read) const
    {
        return !read.is_marked_supplementary_alignment();
    }
};

struct is_good_mapping_quality
{
    using QualityType = AlignedRead::QualityType;
    
    is_good_mapping_quality() = default;
    explicit is_good_mapping_quality(QualityType good_mapping_quality) : good_mapping_quality_ {good_mapping_quality} {}
    
    bool operator()(const AlignedRead& read) const
    {
        return read.get_mapping_quality() >= good_mapping_quality_;
    }
    
private:
    const QualityType good_mapping_quality_;
};

struct has_good_base_fraction
{
    using QualityType = AlignedRead::QualityType;
    
    has_good_base_fraction() = default;
    explicit has_good_base_fraction(QualityType good_base_quality, double min_good_base_fraction)
    : good_base_quality_ {good_base_quality}, min_good_base_fraction_ {min_good_base_fraction} {}
    
    bool operator()(const AlignedRead& read) const
    {
        const auto& qualities = read.get_qualities();
        auto num_good_bases = std::count_if(std::cbegin(qualities), std::cend(qualities), [this]
                                            (auto quality) { return quality >= good_base_quality_; });
        auto good_base_fraction = static_cast<double>(num_good_bases) / static_cast<double>(sequence_size(read));
        return good_base_fraction >= min_good_base_fraction_;
    }
    
private:
    const QualityType good_base_quality_;
    double min_good_base_fraction_;
};

struct has_sufficient_good_quality_bases
{
    using QualityType = AlignedRead::QualityType;
    
    has_sufficient_good_quality_bases() = default;
    explicit has_sufficient_good_quality_bases(QualityType good_base_quality, unsigned min_good_bases)
    : good_base_quality_ {good_base_quality}, min_good_bases_ {min_good_bases} {}
    
    bool operator()(const AlignedRead& read) const
    {
        const auto& qualities = read.get_qualities();
        return std::count_if(std::cbegin(qualities), std::cend(qualities), [this]
                             (auto quality) { return quality >= good_base_quality_; }) >= min_good_bases_;
    }
    
private:
    const QualityType good_base_quality_;
    const unsigned min_good_bases_;
};

struct is_mapped
{
    is_mapped() = default;
    
    bool operator()(const AlignedRead& read) const
    {
        return !read.is_marked_unmapped();
    }
};

struct is_not_chimeric
{
    is_not_chimeric() = default;
    
    bool operator()(const AlignedRead& read) const
    {
        return !read.is_chimeric();
    }
};

//struct is_next_segment_mapped
//{
//    is_next_segment_mapped() = default;
//    
//    bool operator()(const AlignedRead& read) const
//    {
//        return !read.has_mate_pair() || read.is_marked_proper_pair();
//    }
//};

struct is_not_marked_duplicate
{
    is_not_marked_duplicate() = default;
    
    bool operator()(const AlignedRead& read) const
    {
        return !read.is_marked_duplicate();
    }
};

struct is_short
{
    using SizeType = AlignedRead::SizeType;
    
    is_short() = default;
    explicit is_short(SizeType max_length) : max_length_ {max_length} {}
    
    bool operator()(const AlignedRead& read) const
    {
        return sequence_size(read) <= max_length_;
    }
    
private:
    const SizeType max_length_;
};

struct is_long
{
    using SizeType = AlignedRead::SizeType;
    
    is_long() = default;
    explicit is_long(SizeType min_length) : min_length_ {min_length} {}
    
    bool operator()(const AlignedRead& read) const
    {
        return read.get_mapping_quality() >= min_length_;
    }
    
private:
    const SizeType min_length_;
};

struct is_not_contaminated
{
    is_not_contaminated() = default;
    
    bool operator()(const AlignedRead& read) const
    {
        return !read.is_chimeric() || sequence_size(read) >= read.get_next_segment().get_inferred_template_length();
    }
};

struct is_not_marked_qc_fail
{
    is_not_marked_qc_fail() = default;
    
    bool operator()(const AlignedRead& read) const
    {
        return !read.is_marked_qc_fail();
    }
};

struct mate_is_mapped
{
    mate_is_mapped() = default;
    
    bool operator()(const AlignedRead& read) const
    {
        return !read.has_mate() || !read.get_next_segment().is_marked_unmapped();
    }
};

// Context-based filters

struct is_not_duplicate
{
    is_not_duplicate() = default;
    
    template <typename Iterator>
    bool operator()(const AlignedRead& read, Iterator first_good_read, Iterator last_good_read) const
    {
        return first_good_read == last_good_read || !IsDuplicate()(read, *last_good_read);
    }
};

} // namespace ReadFilters
} // namespace Octopus

#endif /* defined(__Octopus__read_filters__) */
