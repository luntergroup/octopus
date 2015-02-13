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

#include "genomic_region.h"
#include "cigar_string.h"

using std::uint_fast32_t;
using std::uint_fast8_t;

class AlignedRead
{
public:
    AlignedRead() = delete;
    AlignedRead(GenomicRegion reference_region, std::string sequence,
                std::vector<uint_fast8_t> qualities, CigarString cigar_string,
                uint_fast32_t insert_size, std::string mate_contig_name,
                uint_fast32_t mate_begin, uint_fast8_t mapping_quality);
    
    GenomicRegion get_region() const;
    const std::string& get_contig_name() const;
    uint_fast32_t get_begin() const;
    uint_fast32_t get_end() const;
    const std::string& get_sequence() const;
    const std::vector<uint_fast8_t>& get_qualities() const;
    uint_fast32_t get_sequence_size() const;
    const CigarString& get_cigar_string() const;
    uint_fast32_t get_insert_size() const;
    const std::string& get_mate_contig_name() const;
    uint_fast32_t get_mate_begin() const;
    
private:
    const std::string sequence_;
    const std::vector<uint_fast8_t> qualities_;
    const CigarString cigar_string_;
    const GenomicRegion reference_region_;
    const uint_fast32_t insert_size_;
    const std::string mate_contig_name_;
    const uint_fast32_t mate_begin_;
    const uint_fast8_t mapping_quality_;
};

inline
AlignedRead::AlignedRead(GenomicRegion reference_region, std::string sequence,
                         std::vector<uint_fast8_t> qualities, CigarString cigar_string,
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

inline
GenomicRegion AlignedRead::get_region() const
{
    return reference_region_;
}

inline
const std::string& AlignedRead::get_contig_name() const
{
    return reference_region_.get_contig_name();
}

inline
uint_fast32_t AlignedRead::get_begin() const
{
    return reference_region_.get_begin();
}

inline
uint_fast32_t AlignedRead::get_end() const
{
    return reference_region_.get_end();
}

inline
const std::string& AlignedRead::get_sequence() const
{
    return sequence_;
}

inline
const std::vector<uint_fast8_t>& AlignedRead::get_qualities() const
{
    return qualities_;
}

inline
uint_fast32_t AlignedRead::get_sequence_size() const
{
    return static_cast<uint_fast32_t>(sequence_.size());
}

inline
const CigarString& AlignedRead::get_cigar_string() const
{
    return cigar_string_;
}

inline
uint_fast32_t AlignedRead::get_insert_size() const
{
    return insert_size_;
}

inline
const std::string& AlignedRead::get_mate_contig_name() const
{
    return mate_contig_name_;
}

inline
uint_fast32_t AlignedRead::get_mate_begin() const
{
    return mate_begin_;
}

inline bool operator==(const AlignedRead& lhs, const AlignedRead& rhs)
{
    return lhs.get_region() == rhs.get_region() &&
            lhs.get_sequence() == rhs.get_sequence() &&
            lhs.get_qualities() == rhs.get_qualities() &&
            lhs.get_insert_size() == rhs.get_insert_size();
}
inline bool operator< (const AlignedRead& lhs, const AlignedRead& rhs)
{
    return lhs.get_begin() < rhs.get_begin();
}
inline bool operator!=(const AlignedRead& lhs, const AlignedRead& rhs){return !operator==(lhs, rhs);}
inline bool operator> (const AlignedRead& lhs, const AlignedRead& rhs){return operator<(rhs,lhs);}
inline bool operator<=(const AlignedRead& lhs, const AlignedRead& rhs){return !operator>(lhs,rhs);}
inline bool operator>=(const AlignedRead& lhs, const AlignedRead& rhs){return !operator<(lhs,rhs);}

namespace std {
    template <> struct hash<AlignedRead>
    {
        size_t operator()(const AlignedRead& r) const
        {
            return hash<std::string>()(r.get_sequence());
        }
    };
}

inline
std::ostream& operator<<(std::ostream& os, const AlignedRead& a_read)
{
    os << a_read.get_region() << '\n';
    os << a_read.get_sequence() << '\n';
    //os << a_read.get_qualities() << '\n';
    os << a_read.get_cigar_string() << '\n';
    os << a_read.get_insert_size() << '\n';
    os << a_read.get_mate_contig_name() << '\n';
    os << a_read.get_mate_begin() << '\n';
    return os;
}

#endif /* defined(__Octopus__aligned_read__) */
