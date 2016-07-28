//
//  aligned_read.hpp
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__aligned_read__
#define __Octopus__aligned_read__

#include <string>
#include <cstdint>
#include <cstddef>
#include <vector>
#include <bitset>
#include <algorithm>
#include <iterator>
#include <utility>
#include <functional>
#include <iosfwd>

#include <boost/optional.hpp>

#include "comparable.hpp"
#include "equitable.hpp"
#include "genomic_region.hpp"
#include "mappable.hpp"
#include "cigar_string.hpp"

class AlignedRead : public Comparable<AlignedRead>, public Mappable<AlignedRead>
{
public:
    using RegionType          = GenomicRegion;
    using NucleotideSequence  = std::string;
    using MappingQuality      = std::uint_fast8_t;
    using BaseQuality         = std::uint_fast8_t;
    using BaseQualityVector   = std::vector<BaseQuality>;
    
    class NextSegment : public Equitable<NextSegment>
    {
    public:
        struct Flags;
        
        NextSegment() = default;
        
        template <typename S>
        NextSegment(S&& contig_name, GenomicRegion::Position begin,
                    GenomicRegion::Size inferred_template_length,
                    Flags data);
        
        NextSegment(const NextSegment&)            = default;
        NextSegment& operator=(const NextSegment&) = default;
        NextSegment(NextSegment&&)                 = default;
        NextSegment& operator=(NextSegment&&)      = default;
        
        ~NextSegment() = default;
        
        const GenomicRegion::ContigName& contig_name() const;
        
        GenomicRegion::Position begin() const noexcept;
        
        GenomicRegion::Size inferred_template_length() const noexcept;
        
        bool is_marked_unmapped() const;
        bool is_marked_reverse_mapped() const;
        
        std::size_t get_hash() const;
        
    private:
        using FlagBits = std::bitset<2>;
        
        GenomicRegion::ContigName contig_name_;
        GenomicRegion::Position begin_;
        GenomicRegion::Size inferred_template_length_;
        FlagBits flags_;
        
        FlagBits compress(const Flags& data);
    };
    
    struct Flags;
    
    AlignedRead() = default;
    
    template <typename GenomicRegion_, typename String1_, typename Qualities_, typename CigarString_>
    AlignedRead(GenomicRegion_&& reference_region, String1_&& sequence,
                Qualities_&& qualities, CigarString_&& cigar_string,
                MappingQuality mapping_quality, const Flags& flags);
    
    template <typename GenomicRegion_, typename String1_, typename Qualities_, typename CigarString_,
              typename String2_>
    AlignedRead(GenomicRegion_&& reference_region, String1_&& sequence,
                Qualities_&& qualities, CigarString_&& cigar_string,
                MappingQuality mapping_quality, Flags flags,
                String2_&& next_segment_contig_name, RegionType::Position next_segment_begin,
                RegionType::Size inferred_template_length,
                const NextSegment::Flags& next_segment_flags);
    
    AlignedRead(const AlignedRead& other)            = default;
    AlignedRead& operator=(const AlignedRead& other) = default;
    AlignedRead(AlignedRead&&)                       = default;
    AlignedRead& operator=(AlignedRead&&)            = default;
    
    ~AlignedRead() = default;
    
    const std::string& read_group() const;
    
    const GenomicRegion& mapped_region() const noexcept;
    
    const NucleotideSequence& sequence() const noexcept;
    
    const BaseQualityVector& qualities() const noexcept;
    
    MappingQuality mapping_quality() const noexcept;
    
    const CigarString& cigar_string() const noexcept;
    
    bool has_other_segment() const noexcept;
    
    const NextSegment& next_segment() const;
    
    Flags flags() const;
    
    bool is_marked_all_segments_in_read_aligned() const noexcept;
    bool is_marked_multiple_segment_template() const noexcept;
    bool is_marked_unmapped() const noexcept;
    bool is_marked_reverse_mapped() const noexcept;
    bool is_marked_secondary_alignment() const noexcept;
    bool is_marked_qc_fail() const noexcept;
    bool is_marked_duplicate() const noexcept;
    bool is_marked_supplementary_alignment() const noexcept;
    
    std::size_t get_hash() const;
    
    void cap_qualities(BaseQuality max) noexcept;
    
    void zero_front_qualities(std::size_t num_bases) noexcept;
    void zero_back_qualities(std::size_t num_bases) noexcept;
    
    void compress();
    void decompress();
    
private:
    static constexpr std::size_t Num_flags = 9;
    
    using FlagBits = std::bitset<Num_flags>;
    
    static constexpr std::size_t compression_flag_ {Num_flags - 1};
    
    // should be ordered by sizeof
    GenomicRegion region_;
    NucleotideSequence sequence_;
    BaseQualityVector qualities_;
    CigarString cigar_string_;
    std::string read_group_;
    boost::optional<NextSegment> next_segment_;
    mutable std::size_t hash_ = 0; // 0 is reserved so can be lazy evaluated
    FlagBits flags_;
    MappingQuality mapping_quality_;
    
    FlagBits compress(const Flags& data);
    
    bool is_compressed() const noexcept;
    void set_compressed() noexcept;
    void set_uncompressed() noexcept;
    
    std::size_t make_hash() const;
};

struct AlignedRead::NextSegment::Flags
{
    bool unmapped;
    bool reverse_mapped;
};

struct AlignedRead::Flags
{
    bool multiple_segment_template;
    bool all_segments_in_read_aligned;
    bool unmapped;
    bool reverse_mapped;
    bool first_template_segment;
    bool last_template_segmenet;
    bool secondary_alignment;
    bool qc_fail;
    bool duplicate;
    bool supplementary_alignment;
};

template <typename GenomicRegion_, typename String1_, typename Qualities_, typename CigarString_>
AlignedRead::AlignedRead(GenomicRegion_&& reference_region, String1_&& sequence,
                         Qualities_&& qualities, CigarString_&& cigar_string,
                         MappingQuality mapping_quality, const Flags& flags)
:
region_ {std::forward<GenomicRegion_>(reference_region)},
sequence_ {std::forward<String1_>(sequence)},
qualities_ {std::forward<Qualities_>(qualities)},
cigar_string_ {std::forward<CigarString_>(cigar_string)},
read_group_ {},
next_segment_ {},
flags_ {compress(flags)},
mapping_quality_ {mapping_quality}
{}

template <typename GenomicRegion_, typename String1_, typename Qualities_, typename CigarString_,
          typename String2_>
AlignedRead::AlignedRead(GenomicRegion_&& reference_region, String1_&& sequence,
                         Qualities_&& qualities, CigarString_&& cigar_string,
                         MappingQuality mapping_quality, Flags flags,
                         String2_&& next_segment_contig_name,
                         RegionType::Position next_segment_begin,
                         RegionType::Size inferred_template_length,
                         const NextSegment::Flags& next_segment_flags)
:
region_ {std::forward<GenomicRegion_>(reference_region)},
sequence_ {std::forward<String1_>(sequence)},
qualities_ {std::forward<Qualities_>(qualities)},
cigar_string_ {std::forward<CigarString_>(cigar_string)},
read_group_ {},
next_segment_ {NextSegment {std::forward<String2_>(next_segment_contig_name), next_segment_begin,
    inferred_template_length, next_segment_flags}},
flags_ {compress(flags)},
mapping_quality_ {mapping_quality}
{}

template <typename String_>
AlignedRead::NextSegment::NextSegment(String_&& contig_name, GenomicRegion::Position begin,
                                      GenomicRegion::Size inferred_template_length, Flags data)
:
contig_name_ {std::forward<String_>(contig_name)},
begin_ {begin},
inferred_template_length_ {inferred_template_length},
flags_ {compress(data)}
{}

// Non-member methods

AlignedRead::NucleotideSequence::size_type sequence_size(const AlignedRead& read) noexcept;

bool is_empty_sequence(const AlignedRead& read) noexcept;

CigarString splice_cigar(const AlignedRead& read, const GenomicRegion& region);

ContigRegion::Size count_overlapped_bases(const AlignedRead& read, const GenomicRegion& region);

bool is_soft_clipped(const AlignedRead& read);

std::pair<CigarOperation::Size, CigarOperation::Size> get_soft_clipped_sizes(const AlignedRead& read);

AlignedRead splice(const AlignedRead& read, const GenomicRegion& region);

bool operator==(const AlignedRead& lhs, const AlignedRead& rhs);
bool operator<(const AlignedRead& lhs, const AlignedRead& rhs);

struct IsDuplicate
{
    bool operator()(const AlignedRead& lhs, const AlignedRead& rhs) const;
};

bool operator==(const AlignedRead::NextSegment& lhs, const AlignedRead::NextSegment& rhs);

namespace std {
    template <> struct hash<AlignedRead>
    {
        size_t operator()(const AlignedRead& read) const
        {
            return read.get_hash();
        }
    };
    
    template <> struct hash<reference_wrapper<const AlignedRead>>
    {
        size_t operator()(const reference_wrapper<const AlignedRead> read) const
        {
            return hash<AlignedRead>()(read);
        }
    };
} // namespace std

namespace boost
{
    template <> struct hash<AlignedRead> : std::unary_function<AlignedRead, std::size_t>
    {
        std::size_t operator()(const AlignedRead& read) const
        {
            return std::hash<AlignedRead>()(read);
        }
    };
} // namespace boost

std::ostream& operator<<(std::ostream& os, const AlignedRead::BaseQualityVector& qualities);
std::ostream& operator<<(std::ostream& os, const AlignedRead& read);

#endif /* defined(__Octopus__aligned_read__) */
