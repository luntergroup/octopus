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
#include <algorithm> // std::transform, std::swap
#include <memory>    // std::unique_ptr, std::make_unique
#include <iterator>  // std::begin etc

#include "genomic_region.h"
#include "cigar_string.h"
#include "comparable.h"

class AlignedRead : public Comparable<AlignedRead>
{
public:
    using SizeType    = GenomicRegion::SizeType;
    using StringType  = std::string;
    using QualityType = std::uint_fast8_t;
    using Qualities   = std::vector<QualityType>;
    
    class MatePair : public Equitable<MatePair>
    {
    public:
        struct SupplementaryData
        {
            bool is_marked_unmapped;
            bool is_marked_reverse_mapped;
        };
        
        MatePair() = default;
        template <typename String_>
        MatePair(String_&& contig_name, SizeType begin, SizeType insert_size,
                 SupplementaryData data);
        ~MatePair() = default;
        
        MatePair(const MatePair&)            = default;
        MatePair& operator=(const MatePair&) = default;
        MatePair(MatePair&&)                 = default;
        MatePair& operator=(MatePair&&)      = default;
        
        const std::string& get_contig_name() const;
        SizeType get_begin() const noexcept;
        SizeType get_insert_size() const noexcept;
        
        bool is_marked_unmapped() const;
        bool is_marked_reverse_mapped() const;
        
    private:
        using Flags = std::vector<bool>;
        
        std::string the_contig_name_;
        SizeType begin_;
        SizeType the_insert_size_;
        Flags flags_;
        
        Flags get_flags_(const SupplementaryData& data);
    };
    
    struct SupplementaryData
    {
        bool is_marked_duplicate;
        bool is_marked_unmapped;
        bool is_marked_reverse_mapped;
        bool is_marked_paired;
        bool is_marked_proper_pair;
        bool is_marked_secondary_alignment;
        bool is_marked_qc_fail;
    };
    
    AlignedRead() = delete;
    
    // For reads without a mate pair
    template <typename GenomicRegion_, typename String1_, typename Qualities_, typename CigarString_>
    explicit AlignedRead(GenomicRegion_&& reference_region, String1_&& sequence,
                         Qualities_&& qualities, CigarString_&& cigar_string,
                         QualityType mapping_quality, SupplementaryData data);
    
    // For reads with a mate pair
    template <typename GenomicRegion_, typename String1_, typename Qualities_, typename CigarString_,
              typename String2_>
    explicit AlignedRead(GenomicRegion_&& reference_region, String1_&& sequence,
                         Qualities_&& qualities, CigarString_&& cigar_string,
                         SizeType insert_size, String2_&& mate_contig_name,
                         SizeType mate_begin, QualityType mapping_quality,
                         SupplementaryData data,
                         MatePair::SupplementaryData mate_pair_data);
    
    AlignedRead(const AlignedRead& other);
    AlignedRead& operator=(const AlignedRead& other);
    AlignedRead(AlignedRead&&)            = default;
    AlignedRead& operator=(AlignedRead&&) = default;
    friend void swap(AlignedRead& lhs, AlignedRead& rhs) noexcept;
    
    // Things that should never change (in theory)
    const GenomicRegion& get_region() const;
    const std::string& get_contig_name() const;
    SizeType get_begin() const noexcept;
    SizeType get_end() const noexcept;
    const StringType& get_sequence() const;
    const Qualities& get_qualities() const;
    QualityType get_mapping_quality() const;
    SizeType get_sequence_size() const;
    const CigarString& get_cigar_string() const;
    
    void set_qualities(Qualities&& new_qualities) noexcept;
    void zero_front_qualities(SizeType num_bases) noexcept;
    void zero_back_qualities(SizeType num_bases) noexcept;
    
    // Mate-pair stuff
    bool has_mate_pair() const noexcept;
    const std::unique_ptr<MatePair>& get_mate_pair() const;
    
    // Things that could change / dependent on Samtools specification
    bool is_marked_duplicate() const;
    bool is_marked_unmapped() const;
    bool is_marked_reverse_mapped() const;
    bool is_marked_paired() const;
    bool is_marked_proper_pair() const;
    bool is_marked_secondary_alignment() const;
    bool is_marked_qc_fail() const;
    
    template <typename CompressionAlgorithm> void compress(const CompressionAlgorithm& c);
    template <typename CompressionAlgorithm> void decompress(const CompressionAlgorithm& c);
    
private:
    using Flags = std::vector<bool>;
    
    GenomicRegion the_reference_region_;
    StringType the_sequence_;
    CigarString the_cigar_string_;
    Qualities the_qualities_;
    std::unique_ptr<MatePair> the_mate_pair_;
    QualityType the_mapping_quality_;
    Flags flags_;
    
    Flags get_flags_(const SupplementaryData& data);
    
    bool is_compressed_() const noexcept;
    void set_compressed_() noexcept;
    void set_uncompressed_() noexcept;
};

template <typename GenomicRegion_, typename String1_, typename Qualities_, typename CigarString_>
inline AlignedRead::AlignedRead(GenomicRegion_&& reference_region, String1_&& sequence,
                                Qualities_&& qualities, CigarString_&& cigar_string,
                                QualityType mapping_quality, SupplementaryData supplementary_data)
:
the_reference_region_ {std::forward<GenomicRegion_>(reference_region)},
the_sequence_ {std::forward<String1_>(sequence)},
the_qualities_ {std::forward<Qualities_>(qualities)},
the_cigar_string_ {std::forward<CigarString_>(cigar_string)},
the_mate_pair_ {nullptr},
flags_ {get_flags_(supplementary_data)}
{}

template <typename GenomicRegion_, typename String1_, typename Qualities_, typename CigarString_,
          typename String2_>
inline AlignedRead::AlignedRead(GenomicRegion_&& reference_region, String1_&& sequence,
                                Qualities_&& qualities, CigarString_&& cigar_string,
                                SizeType insert_size, String2_&& mate_contig_name,
                                SizeType mate_begin, QualityType mapping_quality,
                                SupplementaryData data,
                                MatePair::SupplementaryData mate_pair_data)
:
the_reference_region_ {std::forward<GenomicRegion_>(reference_region)},
the_sequence_ {std::forward<String1_>(sequence)},
the_qualities_ {std::forward<Qualities_>(qualities)},
the_cigar_string_ {std::forward<CigarString_>(cigar_string)},
the_mate_pair_ {std::make_unique<MatePair>(std::forward<String2_>(mate_contig_name),
                                           mate_begin, insert_size, mate_pair_data)},
the_mapping_quality_ {mapping_quality},
flags_ {get_flags_(data)}
{}

template <typename String_>
AlignedRead::MatePair::MatePair(String_&& contig_name, SizeType begin, SizeType insert_size,
                                SupplementaryData data)
:
the_contig_name_ {std::forward<String_>(contig_name)},
begin_ {begin},
the_insert_size_ {insert_size},
flags_ {get_flags_(data)}
{}

inline AlignedRead::AlignedRead(const AlignedRead& other)
:
the_reference_region_ {other.the_reference_region_},
the_sequence_ {other.the_sequence_},
the_qualities_ {other.the_qualities_},
the_cigar_string_ {other.the_cigar_string_},
the_mate_pair_ {(other.the_mate_pair_ != nullptr) ?
    std::make_unique<MatePair>(*other.the_mate_pair_) : nullptr},
the_mapping_quality_ {other.the_mapping_quality_},
flags_ {other.flags_}
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
    std::swap(lhs.the_sequence_, rhs.the_sequence_);
    std::swap(lhs.the_cigar_string_, rhs.the_cigar_string_);
    std::swap(lhs.the_qualities_, rhs.the_qualities_);
    std::swap(lhs.the_mate_pair_, rhs.the_mate_pair_);
    std::swap(lhs.the_mapping_quality_, rhs.the_mapping_quality_);
    std::swap(lhs.flags_, rhs.flags_);
}

inline const std::string& AlignedRead::MatePair::get_contig_name() const
{
    return the_contig_name_;
}

inline AlignedRead::SizeType AlignedRead::MatePair::get_begin() const noexcept
{
    return begin_;
}

inline AlignedRead::SizeType AlignedRead::MatePair::get_insert_size() const noexcept
{
    return the_insert_size_;
}

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

inline const AlignedRead::StringType& AlignedRead::get_sequence() const
{
    return the_sequence_;
}

inline const AlignedRead::Qualities& AlignedRead::get_qualities() const
{
    return the_qualities_;
}

inline void AlignedRead::set_qualities(Qualities&& new_qualities) noexcept
{
    the_qualities_ = std::move(new_qualities);
}

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
    
inline AlignedRead::QualityType AlignedRead::get_mapping_quality() const
{
    return the_mapping_quality_;
}

inline AlignedRead::SizeType AlignedRead::get_sequence_size() const
{
    return static_cast<SizeType>(the_sequence_.size());
}

inline const CigarString& AlignedRead::get_cigar_string() const
{
    return the_cigar_string_;
}

inline bool AlignedRead::has_mate_pair() const noexcept
{
    return the_mate_pair_ != nullptr;
}

inline const std::unique_ptr<AlignedRead::MatePair>& AlignedRead::get_mate_pair() const
{
    if (has_mate_pair()) {
        return the_mate_pair_;
    } else {
        throw std::runtime_error {"Read does not have a mate-pair"};
    }
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

inline AlignedRead::Flags AlignedRead::get_flags_(const SupplementaryData& data)
{
    Flags result(7);
    result[0] = false;
    result[1] = data.is_marked_duplicate;
    result[2] = data.is_marked_unmapped;
    result[3] = data.is_marked_reverse_mapped;
    result[4] = data.is_marked_paired;
    result[5] = data.is_marked_proper_pair;
    result[6] = data.is_marked_secondary_alignment;
    result[7] = data.is_marked_qc_fail;
    return result;
}

inline AlignedRead::MatePair::Flags
AlignedRead::MatePair::get_flags_(const SupplementaryData& data)
{
    Flags result(2);
    result[0] = data.is_marked_unmapped;
    result[1] = data.is_marked_reverse_mapped;
    return result;
}

inline bool AlignedRead::is_marked_duplicate() const
{
    return flags_[1];
}

inline bool AlignedRead::is_marked_unmapped() const
{
    return flags_[2];
}

inline bool AlignedRead::is_marked_reverse_mapped() const
{
    return flags_[3];
}

inline bool AlignedRead::is_marked_paired() const
{
    return flags_[4];
}

inline bool AlignedRead::is_marked_proper_pair() const
{
    return flags_[5];
}

inline bool AlignedRead::is_marked_secondary_alignment() const
{
    return flags_[6];
}

inline bool AlignedRead::is_marked_qc_fail() const
{
    return flags_[7];
}


inline bool AlignedRead::is_compressed_() const noexcept
{
    return flags_[0];
}

inline void AlignedRead::set_compressed_() noexcept
{
    flags_[0] = true;
}

inline void AlignedRead::set_uncompressed_() noexcept
{
    flags_[0] = false;
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

inline bool operator==(const AlignedRead::MatePair& lhs, const AlignedRead::MatePair& rhs)
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
    if (a_read.has_mate_pair()) {
        os << '\n';
        os << a_read.get_mate_pair()->get_contig_name() << '\n';
        os << a_read.get_mate_pair()->get_begin() << '\n';
        os << a_read.get_mate_pair()->get_insert_size();
    }
    return os;
}

#endif /* defined(__Octopus__aligned_read__) */
