//
//  read_filters.cpp
//  Octopus
//
//  Created by Daniel Cooke on 26/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "read_filters.hpp"

namespace Octopus { namespace ReadFilters
{
    bool IsNotSecondaryAlignment::passes(const AlignedRead& read) const noexcept
    {
        return !read.is_marked_secondary_alignment();
    }
    
    bool IsNotSupplementaryAlignment::passes(const AlignedRead& read) const noexcept
    {
        return !read.is_marked_supplementary_alignment();
    }
    
    IsGoodMappingQuality::IsGoodMappingQuality(QualityType good_mapping_quality)
    :
    BasicReadFilter {"IsGoodMappingQuality"},
    good_mapping_quality_ {good_mapping_quality} {}
    
    bool IsGoodMappingQuality::passes(const AlignedRead& read) const noexcept
    {
        return read.get_mapping_quality() >= good_mapping_quality_;
    }
    
    HasSufficientGoodBaseFraction::HasSufficientGoodBaseFraction(QualityType good_base_quality,
                                                                 double min_good_base_fraction)
    :
    BasicReadFilter {"HasSufficientGoodBaseFraction"},
    good_base_quality_ {good_base_quality}, min_good_base_fraction_ {min_good_base_fraction} {}
    
    bool HasSufficientGoodBaseFraction::passes(const AlignedRead& read) const noexcept
    {
        const auto& qualities = read.get_qualities();
        
        auto num_good_bases = std::count_if(std::cbegin(qualities), std::cend(qualities),
                                            [this] (auto quality) {
                                                return quality >= good_base_quality_;
                                            });
        
        auto good_base_fraction = static_cast<double>(num_good_bases)
                    / static_cast<double>(sequence_size(read));
        
        return good_base_fraction >= min_good_base_fraction_;
    }
    
    HasSufficientGoodQualityBases::HasSufficientGoodQualityBases(QualityType good_base_quality,
                                                                 unsigned min_good_bases)
    :
    BasicReadFilter {"HasSufficientGoodQualityBases"},
    good_base_quality_ {good_base_quality}, min_good_bases_ {min_good_bases} {}
    
    bool HasSufficientGoodQualityBases::passes(const AlignedRead& read) const noexcept
    {
        const auto& qualities = read.get_qualities();
        return std::count_if(std::cbegin(qualities), std::cend(qualities), [this]
                             (auto quality) {
                                 return quality >= good_base_quality_;
                             }) >= min_good_bases_;
    }
    
    bool IsMapped::passes(const AlignedRead& read) const noexcept
    {
        return !read.is_marked_unmapped();
    }
    
    bool IsNotChimeric::passes(const AlignedRead& read) const noexcept
    {
        return !read.is_chimeric();
    }
    
    bool IsNextSegmentMapped::passes(const AlignedRead& read) const noexcept
    {
        return !read.has_mate() || !read.get_next_segment().is_marked_unmapped();
    }
    
    bool IsNotMarkedDuplicate::passes(const AlignedRead& read) const noexcept
    {
        return !read.is_marked_duplicate();
    }
    
    IsShort::IsShort(SizeType max_length)
    :
    BasicReadFilter {"IsShort"},
    max_length_ {max_length} {}
    
    bool IsShort::passes(const AlignedRead& read) const noexcept
    {
        return sequence_size(read) <= max_length_;
    }
    
    IsLong::IsLong(SizeType max_length)
    :
    BasicReadFilter {"IsLong"},
    max_length_ {max_length} {}
    
    bool IsLong::passes(const AlignedRead& read) const noexcept
    {
        return sequence_size(read) >= max_length_;
    }
    
    bool IsNotContaminated::passes(const AlignedRead& read) const noexcept
    {
        return !read.is_chimeric() || sequence_size(read) >= read.get_next_segment().get_inferred_template_length();
    }
    
    bool IsNotMarkedQcFail::passes(const AlignedRead& read) const noexcept
    {
        return !read.is_marked_qc_fail();
    }
} // namespace ReadFilters
} // namespace Octopus