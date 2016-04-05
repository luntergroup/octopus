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
#include <string>

#include "aligned_read.hpp"
#include "cigar_string.hpp"

namespace Octopus { namespace ReadFilters {

// Context-free filters

struct is_not_secondary_alignment
{
    bool operator()(const AlignedRead& read) const
    {
        return !read.is_marked_secondary_alignment();
    }
    
    const std::string& name() const noexcept { return name_; }
    
private:
    std::string name_ = "is_not_secondary_alignment";
};

struct is_not_supplementary_alignment
{
    bool operator()(const AlignedRead& read) const
    {
        return !read.is_marked_supplementary_alignment();
    }
    
    const std::string& name() const noexcept { return name_; }
    
private:
    std::string name_ = "is_not_supplementary_alignment";
};

struct is_good_mapping_quality
{
    using QualityType = AlignedRead::QualityType;
    
    is_good_mapping_quality() = default;
    explicit is_good_mapping_quality(QualityType good_mapping_quality)
    : good_mapping_quality_ {good_mapping_quality} {}
    
    bool operator()(const AlignedRead& read) const
    {
        return read.get_mapping_quality() >= good_mapping_quality_;
    }
    
    const std::string& name() const noexcept { return name_; }
    
private:
    const QualityType good_mapping_quality_;
    
    std::string name_ = "is_good_mapping_quality";
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
    
    const std::string& name() const noexcept { return name_; }
    
private:
    const QualityType good_base_quality_;
    double min_good_base_fraction_;
    
    std::string name_ = "has_good_base_fraction";
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
    
    const std::string& name() const noexcept { return name_; }

private:
    const QualityType good_base_quality_;
    const unsigned min_good_bases_;
    
    std::string name_ = "has_sufficient_good_quality_bases";
};

struct is_mapped
{
    bool operator()(const AlignedRead& read) const
    {
        return !read.is_marked_unmapped();
    }
    
    const std::string& name() const noexcept { return name_; }
    
private:
    std::string name_ = "is_mapped";
};

struct is_not_chimeric
{
    bool operator()(const AlignedRead& read) const
    {
        return !read.is_chimeric();
    }
    
    const std::string& name() const noexcept { return name_; }
    
private:
    std::string name_ = "is_not_chimeric";
};

struct is_next_segment_mapped
{
    bool operator()(const AlignedRead& read) const
    {
        return !read.has_mate() || !read.get_next_segment().is_marked_unmapped();
    }
    
    const std::string& name() const noexcept { return name_; }
    
private:
    std::string name_ = "is_next_segment_mapped";
};

struct is_not_marked_duplicate
{
    bool operator()(const AlignedRead& read) const
    {
        return !read.is_marked_duplicate();
    }
    
    const std::string& name() const noexcept { return name_; }
    
private:
    std::string name_ = "is_not_marked_duplicate";
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
    
    const std::string& name() const noexcept { return name_; }
    
private:
    const SizeType max_length_;
    
    std::string name_ = "is_short";
};

struct is_long
{
    using SizeType = AlignedRead::SizeType;
    
    is_long() = default;
    
    explicit is_long(SizeType max_length) : max_length_ {max_length} {}
    
    bool operator()(const AlignedRead& read) const
    {
        return sequence_size(read) >= max_length_;
    }
    
    const std::string& name() const noexcept { return name_; }
    
private:
    const SizeType max_length_;
    
     std::string name_ = "is_long";
};

struct is_not_contaminated
{
    bool operator()(const AlignedRead& read) const
    {
        return !read.is_chimeric() || sequence_size(read) >= read.get_next_segment().get_inferred_template_length();
    }
    
    const std::string& name() const noexcept { return name_; }
    
private:
    std::string name_ = "is_not_contaminated";
};

struct is_not_marked_qc_fail
{
    bool operator()(const AlignedRead& read) const
    {
        return !read.is_marked_qc_fail();
    }
    
    const std::string& name() const noexcept { return name_; }
    
private:
    std::string name_ = "is_not_marked_qc_fail";
};

// Context-based filters

struct filter_duplicates
{
    template <typename ForwardIt>
    ForwardIt operator()(ForwardIt first_read, ForwardIt last_read) const
    {
        return std::unique(first_read, last_read, IsDuplicate {});
    }
};

} // namespace ReadFilters
} // namespace Octopus

#endif /* defined(__Octopus__read_filters__) */
