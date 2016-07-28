//
//  aligned_read.cpp
//  Octopus
//
//  Created by Daniel Cooke on 18/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "aligned_read.hpp"

#include <ostream>

#include <boost/functional/hash.hpp>

#include "compression.hpp"

//
// NextSegment public methods
//

const std::string& AlignedRead::NextSegment::contig_name() const
{
    return contig_name_;
}

GenomicRegion::Position AlignedRead::NextSegment::begin() const noexcept
{
    return begin_;
}

GenomicRegion::Size AlignedRead::NextSegment::inferred_template_length() const noexcept
{
    return inferred_template_length_;
}

bool AlignedRead::NextSegment::is_marked_unmapped() const
{
    return flags_[0];
}

bool AlignedRead::NextSegment::is_marked_reverse_mapped() const
{
    return flags_[1];
}

//
// AlignedRead public methods
//

const GenomicRegion& AlignedRead::mapped_region() const noexcept
{
    return region_;
}

const AlignedRead::NucleotideSequence& AlignedRead::sequence() const noexcept
{
    return sequence_;
}

const AlignedRead::BaseQualityVector& AlignedRead::qualities() const noexcept
{
    return qualities_;
}

AlignedRead::MappingQuality AlignedRead::mapping_quality() const noexcept
{
    return mapping_quality_;
}

const CigarString& AlignedRead::cigar_string() const noexcept
{
    return cigar_string_;
}

bool AlignedRead::has_other_segment() const noexcept
{
    return static_cast<bool>(next_segment_);
}

const AlignedRead::NextSegment& AlignedRead::next_segment() const
{
    if (has_other_segment()) {
        return *next_segment_;
    } else {
        throw std::runtime_error {"AlignedRead: read does not have a next segment"};
    }
}

AlignedRead::Flags AlignedRead::flags() const
{
    Flags flags {};
    
    flags.multiple_segment_template    = is_marked_multiple_segment_template();
    flags.all_segments_in_read_aligned = is_marked_all_segments_in_read_aligned();
    flags.unmapped                     = is_marked_unmapped();
    flags.reverse_mapped               = is_marked_reverse_mapped();
    flags.secondary_alignment          = is_marked_secondary_alignment();
    flags.qc_fail                      = is_marked_qc_fail();
    flags.duplicate                    = is_marked_duplicate();
    flags.supplementary_alignment      = is_marked_supplementary_alignment();
    
    return flags;
}

bool AlignedRead::is_marked_all_segments_in_read_aligned() const noexcept
{
    return flags_[0];
}

bool AlignedRead::is_marked_multiple_segment_template() const noexcept
{
    return flags_[1];
}

bool AlignedRead::is_marked_unmapped() const noexcept
{
    return flags_[2];
}

bool AlignedRead::is_marked_reverse_mapped() const noexcept
{
    return flags_[3];
}

bool AlignedRead::is_marked_secondary_alignment() const noexcept
{
    return flags_[4];
}

bool AlignedRead::is_marked_qc_fail() const noexcept
{
    return flags_[5];
}

bool AlignedRead::is_marked_duplicate() const noexcept
{
    return flags_[6];
}

bool AlignedRead::is_marked_supplementary_alignment() const noexcept
{
    return flags_[7];
}

std::size_t AlignedRead::get_hash() const
{
    if (hash_ == 0) { // 0 is reserved
        hash_ = make_hash(); // lazy evaluation
    }
    return hash_;
}

void AlignedRead::cap_qualities(const BaseQuality max) noexcept
{
    std::transform(std::cbegin(qualities_), std::cend(qualities_), std::begin(qualities_),
                   [max] (const auto q) { return std::min(q, max); });
}

void AlignedRead::zero_front_qualities(const std::size_t num_bases) noexcept
{
    std::fill_n(std::begin(qualities_), std::min(num_bases, sequence_.size()), 0);
}

void AlignedRead::zero_back_qualities(const std::size_t num_bases) noexcept
{
    std::fill_n(std::rbegin(qualities_), std::min(num_bases, sequence_.size()), 0);
}

void AlignedRead::compress()
{
    sequence_ = octopus::compress(sequence_);
    // TODO: can we also compress qualities and cigar string?
    set_compressed();
}

void AlignedRead::decompress()
{
    sequence_ = octopus::decompress(sequence_);
    set_uncompressed();
}

//
// Private methods
//

bool AlignedRead::is_compressed() const noexcept
{
    return flags_[compression_flag_];
}

void AlignedRead::set_compressed() noexcept
{
    flags_[compression_flag_] = true;
}

void AlignedRead::set_uncompressed() noexcept
{
    flags_[compression_flag_] = false;
}

AlignedRead::FlagBits AlignedRead::compress(const Flags& flags)
{
    FlagBits result {};
    
    result[0] = flags.all_segments_in_read_aligned;
    result[1] = flags.multiple_segment_template;
    result[2] = flags.unmapped;
    result[3] = flags.reverse_mapped;
    result[4] = flags.secondary_alignment;
    result[5] = flags.qc_fail;
    result[6] = flags.duplicate;
    result[7] = flags.supplementary_alignment;
    
    result[compression_flag_] = false;
    
    return result;
}

AlignedRead::NextSegment::FlagBits AlignedRead::NextSegment::compress(const Flags& flags)
{
    FlagBits result {};
    result[0] = flags.unmapped;
    result[1] = flags.reverse_mapped;
    return result;
}

std::size_t AlignedRead::make_hash() const
{
    using boost::hash_combine;
    
    std::size_t result {};
    
    hash_combine(result, std::hash<GenomicRegion>()(region_));
    hash_combine(result, std::hash<CigarString>()(cigar_string_));
    hash_combine(result, boost::hash_range(std::cbegin(qualities_), std::cend(qualities_)));
    hash_combine(result, mapping_quality_);
    
    if (result != 0) return result;
    
    hash_combine(result, std::hash<CigarString>()(cigar_string_));
    
    return (result == 0) ? 1 : result; // 0 is reserved
}

// Non-member methods

AlignedRead::NucleotideSequence::size_type sequence_size(const AlignedRead& read) noexcept
{
    return read.sequence().size();
}

bool is_empty_sequence(const AlignedRead& read) noexcept
{
    return read.sequence().empty();
}

CigarString splice_cigar(const AlignedRead& read, const GenomicRegion& region)
{
    if (contains(region, read)) return read.cigar_string();
    
    const auto splice_region = overlapped_region(read, region);
    
    return splice(read.cigar_string(), static_cast<CigarOperation::Size>(begin_distance(read, splice_region)));
}

ContigRegion::Size count_overlapped_bases(const AlignedRead& read, const GenomicRegion& region)
{
    if (contains(region, read)) {
        return static_cast<ContigRegion::Size>(sequence_size(read));
    }
    
    // TODO: not quite right as doesn't account for indels
    return static_cast<ContigRegion::Size>(std::max(GenomicRegion::Distance {0}, overlap_size(read, region)));
}

bool is_soft_clipped(const AlignedRead& read)
{
    return is_soft_clipped(read.cigar_string());
}

std::pair<CigarOperation::Size, CigarOperation::Size> get_soft_clipped_sizes(const AlignedRead& read)
{
    return get_soft_clipped_sizes(read.cigar_string());
}

AlignedRead splice(const AlignedRead& read, const GenomicRegion& region)
{
    using std::cbegin; using std::next;
    
    if (!overlaps(read, region)) {
        throw std::logic_error {"AlignedRead: trying to splice non-overlapped region"};
    }
    
    if (contains(region, read)) return read;
    
    const auto splice_region = overlapped_region(read, region);
    
    const auto reference_offset = static_cast<CigarOperation::Size>(begin_distance(read, splice_region));
    
    const auto uncontained_cigar_splice = splice_reference(read.cigar_string(), reference_offset);
    
    auto contained_cigar_splice   = splice_reference(read.cigar_string(), reference_offset,
                                                     region_size(splice_region));
    
    const auto sequence_offset = sequence_size(uncontained_cigar_splice);
    const auto sequence_length = sequence_size(contained_cigar_splice);
    
    AlignedRead::NucleotideSequence sequence_splice(next(cbegin(read.sequence()), sequence_offset),
                                                    next(cbegin(read.sequence()),
                                                         sequence_offset + sequence_length));
    
    AlignedRead::BaseQualityVector qualities_splice(next(cbegin(read.qualities()), sequence_offset),
                                                    next(cbegin(read.qualities()),
                                                         sequence_offset + sequence_length));
    
    return AlignedRead {
        splice_region,
        std::move(sequence_splice),
        std::move(qualities_splice),
        std::move(contained_cigar_splice),
        read.mapping_quality(),
        read.flags()
    };
}

bool operator==(const AlignedRead& lhs, const AlignedRead& rhs)
{
    return lhs.mapping_quality() == rhs.mapping_quality()
        && lhs.mapped_region()   == rhs.mapped_region()
        && lhs.cigar_string()    == rhs.cigar_string()
        && lhs.sequence()        == rhs.sequence()
        && lhs.qualities()       == rhs.qualities();
}

bool operator<(const AlignedRead& lhs, const AlignedRead& rhs)
{
    if (lhs.mapped_region() == rhs.mapped_region()) {
        if (lhs.mapping_quality() == rhs.mapping_quality()) {
            if (lhs.cigar_string() == rhs.cigar_string()) {
                if (lhs.sequence() == rhs.sequence()) {
                    return lhs.qualities() < rhs.qualities();
                } else {
                    return lhs.sequence() < rhs.sequence();
                }
            } else {
                return lhs.cigar_string() < rhs.cigar_string();
            }
        } else {
            return lhs.mapping_quality() < rhs.mapping_quality();
        }
    } else {
        return lhs.mapped_region() < rhs.mapped_region();
    }
}

bool are_other_segments_duplicates(const AlignedRead &lhs, const AlignedRead &rhs)
{
    if (lhs.has_other_segment() && rhs.has_other_segment()) {
        return lhs.next_segment() == rhs.next_segment();
    }
    return false;
}

bool IsDuplicate::operator()(const AlignedRead &lhs, const AlignedRead &rhs) const
{
    return lhs.mapped_region() == rhs.mapped_region()
        && lhs.cigar_string() == rhs.cigar_string()
        && lhs.flags().reverse_mapped == rhs.flags().reverse_mapped
        && are_other_segments_duplicates(lhs, rhs);
}

bool operator==(const AlignedRead::NextSegment& lhs, const AlignedRead::NextSegment& rhs)
{
    return lhs.contig_name() == rhs.contig_name()
        && lhs.begin() == rhs.begin()
        && lhs.inferred_template_length() == rhs.inferred_template_length();
}

std::ostream& operator<<(std::ostream& os, const AlignedRead::BaseQualityVector& qualities)
{
    std::transform(std::cbegin(qualities), std::cend(qualities),
                   std::ostream_iterator<AlignedRead::BaseQuality>(os),
                   [] (const auto q) { return static_cast<unsigned>(q + 33); }
                   );
    return os;
}

std::ostream& operator<<(std::ostream& os, const AlignedRead& read)
{
    os << read.mapped_region() << ' ';
    os << read.sequence() << ' ';
    os << read.qualities() << ' ';
    os << read.cigar_string() << ' ';
    os << static_cast<unsigned>(read.mapping_quality()) << ' ';
    if (read.has_other_segment()) {
        os << read.next_segment().contig_name() << ' ';
        os << read.next_segment().begin() << ' ';
        os << read.next_segment().inferred_template_length();
    }
    return os;
}
