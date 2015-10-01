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
#include <ostream>
#include <vector>
#include <bitset>
#include <algorithm> // std::transform, std::swap
#include <memory>    // std::unique_ptr, std::make_unique
#include <iterator>  // std::begin etc
#include <utility>   // std::forward
#include <boost/functional/hash.hpp> // boost::hash_combine

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
        struct FlagData
        {
            FlagData() = default;
            
            bool is_marked_unmapped;
            bool is_marked_reverse_mapped;
        };
        
        NextSegment() = default;
        template <typename String_> NextSegment(String_&& contig_name, SizeType begin,
                                                SizeType inferred_template_length, FlagData data);
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
        
    private:
        using Flags = std::bitset<2>;
        
        std::string contig_name_;
        SizeType begin_;
        SizeType inferred_template_length_;
        Flags flags_;
        
        Flags get_flags(const FlagData& data);
    };
    
    struct FlagData
    {
        FlagData() = default;
        
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
                         QualityType mapping_quality, FlagData flags);
    
    template <typename GenomicRegion_, typename String1_, typename Qualities_, typename CigarString_,
              typename String2_>
    explicit AlignedRead(GenomicRegion_&& reference_region, String1_&& sequence,
                         Qualities_&& qualities, CigarString_&& cigar_string,
                         QualityType mapping_quality, FlagData flags,
                         String2_&& next_segment_contig_name, SizeType next_segment_begin,
                         SizeType inferred_template_length, NextSegment::FlagData next_segment_flags);
    
    AlignedRead(const AlignedRead& other);
    AlignedRead& operator=(const AlignedRead& other);
    AlignedRead(AlignedRead&&)            = default;
    AlignedRead& operator=(AlignedRead&&) = default;
    
    friend void swap(AlignedRead& lhs, AlignedRead& rhs) noexcept;
    
    // basic accessors
    const std::string& get_read_group() const;
    const GenomicRegion& get_region() const noexcept;
    const SequenceType& get_sequence() const noexcept;
    const Qualities& get_qualities() const noexcept;
    QualityType get_mapping_quality() const noexcept;
    SizeType get_sequence_size() const noexcept;
    const CigarString& get_cigar_string() const noexcept;
    const std::unique_ptr<NextSegment>& get_next_segment() const;
    FlagData get_flags() const;
    
    // flags
    bool is_chimeric() const noexcept;
    bool is_marked_all_segments_in_read_aligned() const;
    bool is_marked_multiple_read_template() const;
    bool is_marked_unmapped() const;
    bool is_marked_reverse_mapped() const;
    bool is_marked_secondary_alignment() const;
    bool is_marked_qc_fail() const;
    bool is_marked_duplicate() const;
    bool is_marked_supplementary_alignment() const;
    
    // mutables
    void zero_front_qualities(SizeType num_bases) noexcept;
    void zero_back_qualities(SizeType num_bases) noexcept;
    void compress();
    void decompress();
    
private:
    using Flags = std::bitset<8>;
    
    GenomicRegion reference_region_;
    std::string read_group_;
    SequenceType sequence_;
    CigarString cigar_string_;
    Qualities qualities_;
    std::unique_ptr<NextSegment> next_segment_;
    Flags flags_;
    QualityType mapping_quality_;
    bool is_compressed_;
    
    Flags get_flags(const FlagData& data);
    bool is_compressed() const noexcept;
    void set_compressed() noexcept;
    void set_uncompressed() noexcept;
};

template <typename GenomicRegion_, typename String1_, typename Qualities_, typename CigarString_>
inline AlignedRead::AlignedRead(GenomicRegion_&& reference_region, String1_&& sequence,
                                Qualities_&& qualities, CigarString_&& cigar_string,
                                QualityType mapping_quality, FlagData flags)
:
reference_region_ {std::forward<GenomicRegion_>(reference_region)},
read_group_ {},
sequence_ {std::forward<String1_>(sequence)},
qualities_ {std::forward<Qualities_>(qualities)},
cigar_string_ {std::forward<CigarString_>(cigar_string)},
next_segment_ {nullptr},
flags_ {get_flags(flags)},
mapping_quality_ {mapping_quality},
is_compressed_ {false}
{}

template <typename GenomicRegion_, typename String1_, typename Qualities_, typename CigarString_,
          typename String2_>
AlignedRead::AlignedRead(GenomicRegion_&& reference_region, String1_&& sequence,
                         Qualities_&& qualities, CigarString_&& cigar_string,
                         QualityType mapping_quality, FlagData flags,
                         String2_&& next_segment_contig_name, SizeType next_segment_begin,
                         SizeType inferred_template_length, NextSegment::FlagData next_segment_flags)
:
reference_region_ {std::forward<GenomicRegion_>(reference_region)},
read_group_ {},
sequence_ {std::forward<String1_>(sequence)},
qualities_ {std::forward<Qualities_>(qualities)},
cigar_string_ {std::forward<CigarString_>(cigar_string)},
next_segment_ {std::make_unique<NextSegment>(std::forward<String2_>(next_segment_contig_name),
                                                 next_segment_begin, inferred_template_length,
                                                 next_segment_flags)},
flags_ {get_flags(flags)},
mapping_quality_ {mapping_quality},
is_compressed_ {false}
{}

template <typename String_>
AlignedRead::NextSegment::NextSegment(String_&& contig_name, SizeType begin,
                                      SizeType inferred_template_length, FlagData data)
:
contig_name_ {std::forward<String_>(contig_name)},
begin_ {begin},
inferred_template_length_ {inferred_template_length},
flags_ {get_flags(data)}
{}

// Non-member methods

AlignedRead splice(const AlignedRead& read, const GenomicRegion& region);

bool operator==(const AlignedRead& lhs, const AlignedRead& rhs);

bool operator<(const AlignedRead& lhs, const AlignedRead& rhs);

struct IsDuplicate
{
    bool operator()(const AlignedRead& lhs, const AlignedRead& rhs);
};

bool operator==(const AlignedRead::NextSegment& lhs, const AlignedRead::NextSegment& rhs);

namespace std {
    template <> struct hash<AlignedRead>
    {
        size_t operator()(const AlignedRead& r) const
        {
            size_t seed {};
            boost::hash_combine(seed, hash<GenomicRegion>()(r.get_region()));
            boost::hash_combine(seed, hash<CigarString>()(r.get_cigar_string()));
            boost::hash_combine(seed, r.get_mapping_quality());
            return seed;
        }
    };
} // end namespace std

std::ostream& operator<<(std::ostream& os, const AlignedRead::Qualities& qualities);

std::ostream& operator<<(std::ostream& os, const AlignedRead& a_read);

#endif /* defined(__Octopus__aligned_read__) */
