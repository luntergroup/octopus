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
#include <algorithm>

#include "genomic_region.h"
#include "cigar_string.h"
#include "comparable.h"

using std::uint_fast32_t;
using std::uint_fast8_t;

class AlignedRead : Comparable<AlignedRead>
{
public:
    using QualityType = uint_fast8_t;
    using Qualities   = std::vector<QualityType>;
    
    AlignedRead() = delete;
    AlignedRead(GenomicRegion reference_region, std::string sequence,
                Qualities qualities, CigarString cigar_string,
                uint_fast32_t insert_size, std::string mate_contig_name,
                uint_fast32_t mate_begin, uint_fast8_t mapping_quality);
    
    AlignedRead(const AlignedRead&)            = default;
    AlignedRead& operator=(const AlignedRead&) = default;
    AlignedRead(AlignedRead&&)                 = default;
    AlignedRead& operator=(AlignedRead&&)      = default;
    
    const GenomicRegion& get_region() const;
    const std::string& get_contig_name() const;
    uint_fast32_t get_begin() const;
    uint_fast32_t get_end() const;
    const std::string& get_sequence() const;
    const Qualities& get_qualities() const;
    uint_fast8_t get_mapping_quality() const;
    uint_fast32_t get_sequence_size() const;
    const CigarString& get_cigar_string() const;
    uint_fast32_t get_insert_size() const;
    const std::string& get_mate_contig_name() const;
    uint_fast32_t get_mate_begin() const;
    
private:
    // The seemingly arbitary order of these members is for memory alignment optimisation
    GenomicRegion reference_region_;
    std::string mate_contig_name_;
    std::string sequence_;
    CigarString cigar_string_;
    Qualities qualities_;
    uint_fast32_t insert_size_;
    uint_fast32_t mate_begin_;
    uint_fast8_t mapping_quality_;
};

inline AlignedRead::AlignedRead(GenomicRegion reference_region, std::string sequence,
                                Qualities qualities, CigarString cigar_string,
                                uint_fast32_t insert_size, std::string mate_contig_name,
                                uint_fast32_t mate_begin, uint_fast8_t mapping_quality)
:   reference_region_ {std::move(reference_region)},
    sequence_ {std::move(sequence)},
    qualities_ {std::move(qualities)},
    cigar_string_ {std::move(cigar_string)},
    insert_size_ {insert_size},
    mate_contig_name_ {mate_contig_name},
    mate_begin_ {mate_begin},
    mapping_quality_ {mapping_quality}
{}

inline const GenomicRegion& AlignedRead::get_region() const
{
    return reference_region_;
}

inline const std::string& AlignedRead::get_contig_name() const
{
    return reference_region_.get_contig_name();
}

inline uint_fast32_t AlignedRead::get_begin() const
{
    return reference_region_.get_begin();
}

inline uint_fast32_t AlignedRead::get_end() const
{
    return reference_region_.get_end();
}

inline const std::string& AlignedRead::get_sequence() const
{
    return sequence_;
}

inline const AlignedRead::Qualities& AlignedRead::get_qualities() const
{
    return qualities_;
}

inline uint_fast8_t AlignedRead::get_mapping_quality() const
{
    return mapping_quality_;
}

inline uint_fast32_t AlignedRead::get_sequence_size() const
{
    return static_cast<uint_fast32_t>(sequence_.size());
}

inline const CigarString& AlignedRead::get_cigar_string() const
{
    return cigar_string_;
}

inline uint_fast32_t AlignedRead::get_insert_size() const
{
    return insert_size_;
}

inline const std::string& AlignedRead::get_mate_contig_name() const
{
    return mate_contig_name_;
}

inline uint_fast32_t AlignedRead::get_mate_begin() const
{
    return mate_begin_;
}

inline bool operator==(const AlignedRead& lhs, const AlignedRead& rhs)
{
    // The order of these comparisons should ensure optimal lazy evaluation
    return lhs.get_mapping_quality() == rhs.get_mapping_quality() &&
            lhs.get_region() == rhs.get_region() &&
            lhs.get_cigar_string() == rhs.get_cigar_string() &&
            lhs.get_sequence() == rhs.get_sequence() &&
            lhs.get_qualities() == rhs.get_qualities();
}

inline bool operator< (const AlignedRead& lhs, const AlignedRead& rhs)
{
    return lhs.get_begin() < rhs.get_begin();
}

namespace std {
    template <> struct hash<AlignedRead>
    {
        size_t operator()(const AlignedRead& r) const
        {
            return hash<std::string>()(r.get_sequence());
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
    os << static_cast<unsigned>(a_read.get_mapping_quality()) << '\n';
    os << a_read.get_insert_size() << '\n';
    os << a_read.get_mate_contig_name() << '\n';
    os << a_read.get_mate_begin();
    return os;
}

#endif /* defined(__Octopus__aligned_read__) */
