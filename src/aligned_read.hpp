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
#include <ostream>
#include <vector>
#include <bitset>
#include <algorithm>
#include <memory>
#include <iterator>
#include <utility>

#include "genomic_region.hpp"
#include "cigar_string.hpp"
#include "comparable.hpp"
#include "equitable.hpp"
#include "mappable.hpp"

class AlignedRead : public Comparable<AlignedRead>, public Mappable<AlignedRead>
{
public:
    using SizeType     = GenomicRegion::SizeType;
    using SequenceType = std::string;
    using QualityType  = std::uint_fast8_t;
    using Qualities    = std::vector<QualityType>;
    
    class NextSegment : public Equitable<NextSegment>
    {
    public:
        struct Flags
        {
            Flags() = default;
            
            bool is_marked_unmapped;
            bool is_marked_reverse_mapped;
        };
        
        NextSegment() = default;
        template <typename String_> NextSegment(String_&& contig_name, SizeType begin,
                                                SizeType inferred_template_length, Flags data);
        ~NextSegment() = default;
        
        NextSegment(const NextSegment&)            = default;
        NextSegment& operator=(const NextSegment&) = default;
        NextSegment(NextSegment&&)                 = default;
        NextSegment& operator=(NextSegment&&)      = default;
        
        const std::string& get_contig_name() const;
        SizeType get_begin() const noexcept;
        SizeType get_inferred_template_length() const noexcept;
        GenomicRegion get_inferred_region() const;
        bool is_marked_unmapped() const;
        bool is_marked_reverse_mapped() const;
        
        size_t get_hash() const;
        
    private:
        using FlagBits = std::bitset<2>;
        
        std::string contig_name_;
        SizeType begin_;
        SizeType inferred_template_length_;
        FlagBits flags_;
        
        FlagBits compress_flags(const Flags& data);
    };
    
    struct Flags
    {
        Flags() = default;
        
        bool is_marked_multiple_read_template;
        bool is_marked_all_segments_in_read_aligned;
        bool is_marked_unmapped;
        bool is_marked_reverse_mapped;
        bool is_marked_first_template_segment;
        bool is_marked_last_template_segmenet;
        bool is_marked_secondary_alignment;
        bool is_marked_qc_fail;
        bool is_marked_duplicate;
        bool is_marked_supplementary_alignment;
    };
    
    AlignedRead()  = default;
    ~AlignedRead() = default;
    
    template <typename GenomicRegion_, typename String1_, typename Qualities_, typename CigarString_>
    explicit AlignedRead(GenomicRegion_&& reference_region, String1_&& sequence,
                         Qualities_&& qualities, CigarString_&& cigar_string,
                         QualityType mapping_quality, Flags flags);
    
    template <typename GenomicRegion_, typename String1_, typename Qualities_, typename CigarString_,
              typename String2_>
    explicit AlignedRead(GenomicRegion_&& reference_region, String1_&& sequence,
                         Qualities_&& qualities, CigarString_&& cigar_string,
                         QualityType mapping_quality, Flags flags,
                         String2_&& next_segment_contig_name, SizeType next_segment_begin,
                         SizeType inferred_template_length, NextSegment::Flags next_segment_flags);
    
    AlignedRead(const AlignedRead& other);
    AlignedRead& operator=(const AlignedRead& other);
    AlignedRead(AlignedRead&&)            = default;
    AlignedRead& operator=(AlignedRead&&) = default;
    
    friend void swap(AlignedRead& lhs, AlignedRead& rhs) noexcept;
    
    const std::string& get_read_group() const;
    const GenomicRegion& get_region() const noexcept;
    const SequenceType& get_sequence() const noexcept;
    const Qualities& get_qualities() const noexcept;
    QualityType get_mapping_quality() const noexcept;
    SizeType get_sequence_size() const noexcept;
    const CigarString& get_cigar_string() const noexcept;
    bool has_mate() const noexcept;
    const NextSegment& get_next_segment() const;
    Flags get_flags() const;
    
    bool is_chimeric() const noexcept;
    bool is_marked_all_segments_in_read_aligned() const;
    bool is_marked_multiple_read_template() const;
    bool is_marked_unmapped() const;
    bool is_marked_reverse_mapped() const;
    bool is_marked_secondary_alignment() const;
    bool is_marked_qc_fail() const;
    bool is_marked_duplicate() const;
    bool is_marked_supplementary_alignment() const;
    
    size_t get_hash() const;
    
    // mutables
    void zero_front_qualities(SizeType num_bases) noexcept;
    void zero_back_qualities(SizeType num_bases) noexcept;
    void compress();
    void decompress();
    
private:
    using FlagBits = std::bitset<9>;
    
    static constexpr size_t compression_flag_ {8};
    
    // should be ordered by sizeof
    GenomicRegion region_;
    SequenceType sequence_;
    Qualities qualities_;
    CigarString cigar_string_;
    std::string read_group_;
    std::unique_ptr<NextSegment> next_segment_;
    mutable size_t hash_ = 0; // 0 is reserved so can be lazy evaluated
    FlagBits flags_;
    QualityType mapping_quality_;
    
    FlagBits compress_flags(const Flags& data);
    
    bool is_compressed() const noexcept;
    void set_compressed() noexcept;
    void set_uncompressed() noexcept;
    
    size_t make_hash() const;
};

template <typename GenomicRegion_, typename String1_, typename Qualities_, typename CigarString_>
inline AlignedRead::AlignedRead(GenomicRegion_&& reference_region, String1_&& sequence,
                                Qualities_&& qualities, CigarString_&& cigar_string,
                                QualityType mapping_quality, Flags flags)
:
region_ {std::forward<GenomicRegion_>(reference_region)},
sequence_ {std::forward<String1_>(sequence)},
qualities_ {std::forward<Qualities_>(qualities)},
cigar_string_ {std::forward<CigarString_>(cigar_string)},
read_group_ {},
next_segment_ {nullptr},
flags_ {compress_flags(flags)},
mapping_quality_ {mapping_quality}
{}

template <typename GenomicRegion_, typename String1_, typename Qualities_, typename CigarString_,
          typename String2_>
AlignedRead::AlignedRead(GenomicRegion_&& reference_region, String1_&& sequence,
                         Qualities_&& qualities, CigarString_&& cigar_string,
                         QualityType mapping_quality, Flags flags,
                         String2_&& next_segment_contig_name, SizeType next_segment_begin,
                         SizeType inferred_template_length, NextSegment::Flags next_segment_flags)
:
region_ {std::forward<GenomicRegion_>(reference_region)},
sequence_ {std::forward<String1_>(sequence)},
qualities_ {std::forward<Qualities_>(qualities)},
cigar_string_ {std::forward<CigarString_>(cigar_string)},
read_group_ {},
next_segment_ {std::make_unique<NextSegment>(std::forward<String2_>(next_segment_contig_name),
                                                 next_segment_begin, inferred_template_length,
                                                 next_segment_flags)},
flags_ {compress_flags(flags)},
mapping_quality_ {mapping_quality}
{}

template <typename String_>
AlignedRead::NextSegment::NextSegment(String_&& contig_name, SizeType begin,
                                      SizeType inferred_template_length, Flags data)
:
contig_name_ {std::forward<String_>(contig_name)},
begin_ {begin},
inferred_template_length_ {inferred_template_length},
flags_ {compress_flags(data)}
{}

// Non-member methods

AlignedRead::SizeType num_overlapped_bases(const AlignedRead& read, const GenomicRegion& region);

bool is_soft_clipped(const AlignedRead& read);

std::pair<AlignedRead::SizeType, AlignedRead::SizeType> get_soft_clipped_sizes(const AlignedRead& read);

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
} // namespace std

namespace boost
{
    template <> struct hash<AlignedRead> : std::unary_function<AlignedRead, size_t>
    {
        size_t operator()(const AlignedRead& read) const
        {
            return std::hash<AlignedRead>()(read);
        }
    };
} // namespace boost

std::ostream& operator<<(std::ostream& os, const AlignedRead::Qualities& qualities);
std::ostream& operator<<(std::ostream& os, const AlignedRead& read);

#endif /* defined(__Octopus__aligned_read__) */
