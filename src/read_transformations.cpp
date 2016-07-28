//
//  read_transformations.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "read_transformations.hpp"

#include <algorithm>
#include <numeric>
#include <tuple>

namespace octopus
{
namespace ReadTransforms
{
    CapBaseQualities::CapBaseQualities(BaseQuality max) : max_ {max} {}
    
    void CapBaseQualities::operator()(AlignedRead& read) const noexcept
    {
        read.cap_qualities(max_);
    }
    
    void MaskOverlappedSegment::operator()(AlignedRead& read) const noexcept
    {
        // Only reads in the forward direction are masked to prevent double masking
        if (read.has_other_segment() && contig_name(read) == read.next_segment().contig_name()
            && !read.next_segment().is_marked_unmapped() && !read.is_marked_reverse_mapped()) {
            const auto next_segment_begin = read.next_segment().begin();
            
            if (next_segment_begin < mapped_end(read)) {
                const auto overlapped_size = mapped_end(read) - next_segment_begin;
                read.zero_back_qualities(overlapped_size);
            }
        }
    }
    
    void MaskAdapters::operator()(AlignedRead& read) const noexcept
    {
        if (read.has_other_segment() && contig_name(read) == read.next_segment().contig_name()) {
            const auto insert_size = read.next_segment().inferred_template_length();
            const auto read_size   = sequence_size(read);
            
            if (insert_size < read_size) {
                const auto num_adapter_bases = read_size - insert_size;
                
                if (read.is_marked_reverse_mapped()) {
                    read.zero_front_qualities(num_adapter_bases);
                } else {
                    read.zero_back_qualities(num_adapter_bases);
                }
            }
        }
    }
    
    MaskTail::MaskTail(Length num_bases) : num_bases_ {num_bases} {};
    
    void MaskTail::operator()(AlignedRead& read) const noexcept
    {
        if (read.is_marked_reverse_mapped()) {
            read.zero_front_qualities(num_bases_);
        } else {
            read.zero_back_qualities(num_bases_);
        }
    }
    
    void MaskSoftClipped::operator()(AlignedRead& read) const noexcept
    {
        if (is_soft_clipped(read.cigar_string())) {
            const auto soft_clipped_sizes = get_soft_clipped_sizes(read.cigar_string());
            read.zero_front_qualities(soft_clipped_sizes.first);
            read.zero_back_qualities(soft_clipped_sizes.second);
        }
    }
    
    MaskSoftClippedBoundries::MaskSoftClippedBoundries(Length num_bases) : num_bases_ {num_bases} {};
    
    void MaskSoftClippedBoundries::operator()(AlignedRead& read) const noexcept
    {
        if (is_soft_clipped(read)) {
            Length num_front_bases, num_back_bases;
            std::tie(num_front_bases, num_back_bases) = get_soft_clipped_sizes(read);
            
            if (num_front_bases > 0) {
                read.zero_front_qualities(num_front_bases + num_bases_);
            }
            
            if (num_back_bases > 0) {
                read.zero_back_qualities(num_back_bases + num_bases_);
            }
        }
    }
    
    void QualityAdjustedSoftClippedMasker::operator()(AlignedRead& read) const noexcept
    {
        if (is_soft_clipped(read)) {
            using std::cbegin; using std::crbegin; using std::next;
            using std::accumulate; using std::min_element; using std::min;
            
            using S = AlignedRead::NucleotideSequence::size_type;
            
            S num_front_bases, num_back_bases;
            
            std::tie(num_front_bases, num_back_bases) = get_soft_clipped_sizes(read);
            
            const auto& qualities = read.qualities();
            
            using Q = AlignedRead::BaseQuality;
            
            if (num_front_bases > 0) {
                const auto sum = accumulate(cbegin(qualities), next(cbegin(qualities), num_front_bases), 0.0);
                
                const auto mean = static_cast<Q>(sum / num_front_bases);
                
                const auto min_quality = *min_element(cbegin(qualities), next(cbegin(qualities)));
                
                const auto mask_size = num_front_bases + min(static_cast<S>(mean - min_quality), num_front_bases);
                
                read.zero_front_qualities(mask_size);
            }
            
            if (num_back_bases > 0) {
                const auto sum = accumulate(crbegin(qualities), next(crbegin(qualities), num_back_bases), 0.0);
                
                const auto mean = static_cast<Q>(sum / num_back_bases);
                
                const auto min_quality = *min_element(cbegin(qualities), next(cbegin(qualities)));
                
                const auto mask_size = num_back_bases + min(static_cast<S>(mean - min_quality), num_back_bases);
                
                read.zero_back_qualities(mask_size);
            }
        }
    }
} // namespace ReadTransforms
} // namespace octopus
