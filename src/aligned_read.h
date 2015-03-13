//
//  aligned_read.h
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

#include "genomic_region.h"
#include "cigar_string.h"
#include "comparable.h"

class AlignedRead : public Comparable<AlignedRead>
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
        bool is_marked_unmapped() const;
        bool is_marked_reverse_mapped() const;
        
    private:
        using Flags = std::bitset<2>;
        
        std::string the_contig_name_;
        SizeType begin_;
        SizeType the_inferred_template_length_;
        Flags the_flags_;
        
        Flags get_flags(const FlagData& data);
    };
    
    struct FlagData
    {
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
    
    const std::string& get_read_group() const;
    const GenomicRegion& get_region() const;
    const std::string& get_contig_name() const;
    SizeType get_begin() const noexcept;
    SizeType get_end() const noexcept;
    const SequenceType& get_sequence() const;
    const Qualities& get_qualities() const;
    QualityType get_mapping_quality() const noexcept;
    SizeType get_sequence_size() const noexcept;
    const CigarString& get_cigar_string() const;
    const std::unique_ptr<NextSegment>& get_next_segment() const;
    
    //void set_qualities(Qualities&& new_qualities) noexcept;
    void zero_front_qualities(SizeType num_bases) noexcept;
    void zero_back_qualities(SizeType num_bases) noexcept;
    
    bool is_chimeric() const noexcept;
    
    bool is_marked_all_segments_in_read_aligned() const;
    bool is_marked_multiple_read_template() const;
    bool is_marked_unmapped() const;
    bool is_marked_reverse_mapped() const;
    bool is_marked_secondary_alignment() const;
    bool is_marked_qc_fail() const;
    bool is_marked_duplicate() const;
    bool is_marked_supplementary_alignment() const;
    
    template <typename CompressionAlgorithm> void compress(const CompressionAlgorithm& c);
    template <typename CompressionAlgorithm> void decompress(const CompressionAlgorithm& c);
    
private:
    using Flags = std::bitset<8>;
    
    GenomicRegion the_reference_region_;
    std::string the_read_group_;
    SequenceType the_sequence_;
    CigarString the_cigar_string_;
    Qualities the_qualities_;
    std::unique_ptr<NextSegment> the_next_segment_;
    Flags the_flags_;
    QualityType the_mapping_quality_;
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
the_reference_region_ {std::forward<GenomicRegion_>(reference_region)},
the_read_group_ {},
the_sequence_ {std::forward<String1_>(sequence)},
the_qualities_ {std::forward<Qualities_>(qualities)},
the_cigar_string_ {std::forward<CigarString_>(cigar_string)},
the_next_segment_ {nullptr},
the_flags_ {get_flags(flags)},
the_mapping_quality_ {mapping_quality},
is_compressed_ {false}
{}

template <typename GenomicRegion_, typename String1_, typename Qualities_, typename CigarString_,
          typename String2_>
inline
AlignedRead::AlignedRead(GenomicRegion_&& reference_region, String1_&& sequence,
                         Qualities_&& qualities, CigarString_&& cigar_string,
                         QualityType mapping_quality, FlagData flags,
                         String2_&& next_segment_contig_name, SizeType next_segment_begin,
                         SizeType inferred_template_length, NextSegment::FlagData next_segment_flags)
:
the_reference_region_ {std::forward<GenomicRegion_>(reference_region)},
the_read_group_ {},
the_sequence_ {std::forward<String1_>(sequence)},
the_qualities_ {std::forward<Qualities_>(qualities)},
the_cigar_string_ {std::forward<CigarString_>(cigar_string)},
the_next_segment_ {std::make_unique<NextSegment>(std::forward<String2_>(next_segment_contig_name),
                                                 next_segment_begin, inferred_template_length,
                                                 next_segment_flags)},
the_flags_ {get_flags(flags)},
the_mapping_quality_ {mapping_quality},
is_compressed_ {false}
{}

template <typename String_>
AlignedRead::NextSegment::NextSegment(String_&& contig_name, SizeType begin,
                                      SizeType inferred_template_length, FlagData data)
:
the_contig_name_ {std::forward<String_>(contig_name)},
begin_ {begin},
the_inferred_template_length_ {inferred_template_length},
the_flags_ {get_flags(data)}
{}

inline AlignedRead::AlignedRead(const AlignedRead& other)
:
the_reference_region_ {other.the_reference_region_},
the_read_group_ {other.the_read_group_},
the_sequence_ {other.the_sequence_},
the_qualities_ {other.the_qualities_},
the_cigar_string_ {other.the_cigar_string_},
the_next_segment_ {((other.the_next_segment_ != nullptr) ?
        std::make_unique<NextSegment>(*other.the_next_segment_) : nullptr) },
the_flags_ {other.the_flags_},
the_mapping_quality_ {other.the_mapping_quality_},
is_compressed_ {other.is_compressed()}
{}

inline AlignedRead& AlignedRead::operator=(const AlignedRead& other)
{
    AlignedRead temp {other};
    swap(*this, temp);
    return *this;
}

inline void swap(AlignedRead& lhs, AlignedRead& rhs) noexcept
{
    std::swap(lhs.the_reference_region_, rhs.the_reference_region_);
    std::swap(lhs.the_read_group_, rhs.the_read_group_);
    std::swap(lhs.the_sequence_, rhs.the_sequence_);
    std::swap(lhs.the_cigar_string_, rhs.the_cigar_string_);
    std::swap(lhs.the_qualities_, rhs.the_qualities_);
    std::swap(lhs.the_next_segment_, rhs.the_next_segment_);
    std::swap(lhs.the_flags_, rhs.the_flags_);
    std::swap(lhs.the_mapping_quality_, rhs.the_mapping_quality_);
    std::swap(lhs.is_compressed_, rhs.is_compressed_);
}

//
// NextSegment public methods
//

inline const std::string& AlignedRead::NextSegment::get_contig_name() const
{
    return the_contig_name_;
}

inline AlignedRead::SizeType AlignedRead::NextSegment::get_begin() const noexcept
{
    return begin_;
}

inline AlignedRead::SizeType AlignedRead::NextSegment::get_inferred_template_length() const noexcept
{
    return the_inferred_template_length_;
}

inline bool AlignedRead::NextSegment::is_marked_unmapped() const
{
    return the_flags_[0];
}

inline bool AlignedRead::NextSegment::is_marked_reverse_mapped() const
{
    return the_flags_[1];
}

//
// AlignedRead public methods
//

inline const GenomicRegion& AlignedRead::get_region() const
{
    return the_reference_region_;
}

inline const std::string& AlignedRead::get_contig_name() const
{
    return the_reference_region_.get_contig_name();
}

inline AlignedRead::SizeType AlignedRead::get_begin() const noexcept
{
    return the_reference_region_.get_begin();
}

inline AlignedRead::SizeType AlignedRead::get_end() const noexcept
{
    return the_reference_region_.get_end();
}

inline const AlignedRead::SequenceType& AlignedRead::get_sequence() const
{
    return the_sequence_;
}

inline const AlignedRead::Qualities& AlignedRead::get_qualities() const
{
    return the_qualities_;
}

//inline void AlignedRead::set_qualities(Qualities&& new_qualities) noexcept
//{
//    the_qualities_ = std::move(new_qualities);
//}

inline void AlignedRead::zero_front_qualities(SizeType num_bases) noexcept
{
    std::for_each(std::begin(the_qualities_), std::begin(the_qualities_) + num_bases,
                  [] (auto& a_quality) {
        a_quality = 0;
    });
}

inline void AlignedRead::zero_back_qualities(SizeType num_bases) noexcept
{
    std::for_each(std::rbegin(the_qualities_), std::rbegin(the_qualities_) + num_bases,
                  [] (auto& a_quality) {
        a_quality = 0;
    });
}

inline AlignedRead::QualityType AlignedRead::get_mapping_quality() const noexcept
{
    return the_mapping_quality_;
}

inline AlignedRead::SizeType AlignedRead::get_sequence_size() const noexcept
{
    return static_cast<SizeType>(the_sequence_.size());
}

inline const CigarString& AlignedRead::get_cigar_string() const
{
    return the_cigar_string_;
}

inline const std::unique_ptr<AlignedRead::NextSegment>& AlignedRead::get_next_segment() const
{
    if (is_chimeric()) {
        return the_next_segment_;
    } else {
        throw std::runtime_error {"Read does not have a next segment"};
    }
}

inline bool AlignedRead::is_chimeric() const noexcept
{
    return the_next_segment_ != nullptr;
}

inline bool AlignedRead::is_marked_all_segments_in_read_aligned() const
{
    return the_flags_[0];
}

inline bool AlignedRead::is_marked_multiple_read_template() const
{
    return the_flags_[1];
}

inline bool AlignedRead::is_marked_unmapped() const
{
    return the_flags_[2];
}

inline bool AlignedRead::is_marked_reverse_mapped() const
{
    return the_flags_[3];
}

inline bool AlignedRead::is_marked_secondary_alignment() const
{
    return the_flags_[4];
}

inline bool AlignedRead::is_marked_qc_fail() const
{
    return the_flags_[5];
}

inline bool AlignedRead::is_marked_duplicate() const
{
    return the_flags_[6];
}

inline bool AlignedRead::is_marked_supplementary_alignment() const
{
    return the_flags_[7];
}

template <typename CompressionAlgorithm>
void AlignedRead::compress(const CompressionAlgorithm& c)
{
    the_sequence_ = CompressionAlgorithm::compress(the_sequence_);
}

template <typename CompressionAlgorithm>
void AlignedRead::decompress(const CompressionAlgorithm& c)
{
    the_sequence_ = CompressionAlgorithm::decompress(the_sequence_);
}

//
// Private methods
//

inline bool AlignedRead::is_compressed() const noexcept
{
    return is_compressed_;
}

inline void AlignedRead::set_compressed() noexcept
{
    is_compressed_ = true;
}

inline void AlignedRead::set_uncompressed() noexcept
{
    is_compressed_ = false;
}

inline AlignedRead::Flags AlignedRead::get_flags(const FlagData& flags)
{
    Flags result {};
    result[0] = flags.is_marked_all_segments_in_read_aligned;
    result[1] = flags.is_marked_multiple_read_template;
    result[2] = flags.is_marked_unmapped;
    result[3] = flags.is_marked_reverse_mapped;
    result[4] = flags.is_marked_secondary_alignment;
    result[5] = flags.is_marked_qc_fail;
    result[6] = flags.is_marked_duplicate;
    result[7] = flags.is_marked_supplementary_alignment;
    return result;
}

inline AlignedRead::NextSegment::Flags
AlignedRead::NextSegment::get_flags(const FlagData& flags)
{
    Flags result {};
    result[0] = flags.is_marked_unmapped;
    result[1] = flags.is_marked_reverse_mapped;
    return result;
}

inline bool operator==(const AlignedRead& lhs, const AlignedRead& rhs)
{
    return lhs.get_mapping_quality() == rhs.get_mapping_quality() &&
            lhs.get_region() == rhs.get_region() &&
            lhs.get_cigar_string() == rhs.get_cigar_string();
}

inline bool operator<(const AlignedRead& lhs, const AlignedRead& rhs)
{
    return lhs.get_region() < rhs.get_region();
}

inline bool operator==(const AlignedRead::NextSegment& lhs, const AlignedRead::NextSegment& rhs)
{
    return lhs.get_contig_name() == rhs.get_contig_name() && lhs.get_begin() == rhs.get_begin();
}

namespace std {
    template <> struct hash<AlignedRead>
    {
        size_t operator()(const AlignedRead& r) const
        {
            return hash<GenomicRegion>()(r.get_region()); // TODO: improve this hash
        }
    };
}

inline std::ostream& operator<<(std::ostream& os, const AlignedRead::Qualities& qualities)
{
    std::transform(std::cbegin(qualities), std::cend(qualities),
                   std::ostream_iterator<AlignedRead::QualityType>(os, ""),
                   [] (auto q) { return q + 33; });
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const AlignedRead& a_read)
{
    os << a_read.get_region() << '\n';
    os << a_read.get_sequence() << '\n';
    os << a_read.get_qualities() << '\n';
    os << a_read.get_cigar_string() << '\n';
    os << static_cast<unsigned>(a_read.get_mapping_quality());
    if (a_read.is_chimeric()) {
        os << '\n';
        os << a_read.get_next_segment()->get_contig_name() << '\n';
        os << a_read.get_next_segment()->get_begin() << '\n';
        os << a_read.get_next_segment()->get_inferred_template_length() << '\n';
    } else {
        os << '\n';
        os << "no mate";
    }
    return os;
}

#endif /* defined(__Octopus__aligned_read__) */
