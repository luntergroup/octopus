//
//  aligned_read.cpp
//  Octopus
//
//  Created by Daniel Cooke on 18/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "aligned_read.hpp"

#include <boost/functional/hash.hpp>

#include "compression.hpp"

AlignedRead::AlignedRead(const AlignedRead& other)
:
region_ {other.region_},
read_group_ {other.read_group_},
sequence_ {other.sequence_},
qualities_ {other.qualities_},
cigar_string_ {other.cigar_string_},
next_segment_ {((other.next_segment_ != nullptr) ? std::make_unique<NextSegment>(*other.next_segment_) : nullptr) },
flags_ {other.flags_},
hash_ {other.hash_},
mapping_quality_ {other.mapping_quality_}
{}

AlignedRead& AlignedRead::operator=(const AlignedRead& other)
{
    AlignedRead temp {other};
    swap(*this, temp);
    return *this;
}

void swap(AlignedRead& lhs, AlignedRead& rhs) noexcept
{
    using std::swap;
    swap(lhs.region_, rhs.region_);
    swap(lhs.read_group_, rhs.read_group_);
    swap(lhs.sequence_, rhs.sequence_);
    swap(lhs.cigar_string_, rhs.cigar_string_);
    swap(lhs.qualities_, rhs.qualities_);
    swap(lhs.next_segment_, rhs.next_segment_);
    swap(lhs.flags_, rhs.flags_);
    swap(lhs.hash_, rhs.hash_);
    swap(lhs.mapping_quality_, rhs.mapping_quality_);
}

//
// NextSegment public methods
//

const std::string& AlignedRead::NextSegment::contig_name() const
{
    return contig_name_;
}

AlignedRead::SizeType AlignedRead::NextSegment::begin() const noexcept
{
    return begin_;
}

AlignedRead::SizeType AlignedRead::NextSegment::inferred_template_length() const noexcept
{
    return inferred_template_length_;
}

GenomicRegion AlignedRead::NextSegment::inferred_region() const
{
    return GenomicRegion {contig_name_, begin_, inferred_template_length_};
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

const AlignedRead::SequenceType& AlignedRead::sequence() const noexcept
{
    return sequence_;
}

const AlignedRead::Qualities& AlignedRead::qualities() const noexcept
{
    return qualities_;
}

AlignedRead::QualityType AlignedRead::mapping_quality() const noexcept
{
    return mapping_quality_;
}

const CigarString& AlignedRead::cigar_string() const noexcept
{
    return cigar_string_;
}

bool AlignedRead::is_chimeric() const noexcept
{
    return next_segment_ != nullptr;
}

const AlignedRead::NextSegment& AlignedRead::next_segment() const
{
    if (is_chimeric()) {
        return *next_segment_;
    } else {
        throw std::runtime_error {"AlignedRead: read does not have a next segment"};
    }
}

AlignedRead::Flags AlignedRead::flags() const
{
    Flags flags {};
    
    flags.is_marked_multiple_read_template       = is_marked_multiple_read_template();
    flags.is_marked_all_segments_in_read_aligned = is_marked_all_segments_in_read_aligned();
    flags.is_marked_unmapped                     = is_marked_unmapped();
    flags.is_marked_reverse_mapped               = is_marked_reverse_mapped();
    flags.is_marked_secondary_alignment          = is_marked_secondary_alignment();
    flags.is_marked_qc_fail                      = is_marked_qc_fail();
    flags.is_marked_duplicate                    = is_marked_duplicate();
    flags.is_marked_supplementary_alignment      = is_marked_supplementary_alignment();
    
    return flags;
}

bool AlignedRead::is_marked_all_segments_in_read_aligned() const noexcept
{
    return flags_[0];
}

bool AlignedRead::is_marked_multiple_read_template() const noexcept
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

void AlignedRead::zero_front_qualities(const SizeType num_bases) noexcept
{
    std::fill_n(std::begin(qualities_),
                std::min(num_bases, static_cast<SizeType>(sequence_.size())),
                0);
}

void AlignedRead::zero_back_qualities(const SizeType num_bases) noexcept
{
    std::fill_n(std::rbegin(qualities_),
                std::min(num_bases, static_cast<SizeType>(sequence_.size())),
                0);
}

void AlignedRead::compress()
{
    sequence_ = Octopus::compress(sequence_);
    // TODO: can we also compress qualities and cigar string?
    set_compressed();
}

void AlignedRead::decompress()
{
    sequence_ = Octopus::decompress(sequence_);
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

AlignedRead::FlagBits AlignedRead::compress_flags(const Flags& flags)
{
    FlagBits result {};
    
    result[0] = flags.is_marked_all_segments_in_read_aligned;
    result[1] = flags.is_marked_multiple_read_template;
    result[2] = flags.is_marked_unmapped;
    result[3] = flags.is_marked_reverse_mapped;
    result[4] = flags.is_marked_secondary_alignment;
    result[5] = flags.is_marked_qc_fail;
    result[6] = flags.is_marked_duplicate;
    result[7] = flags.is_marked_supplementary_alignment;
    
    result[compression_flag_] = false;
    
    return result;
}

AlignedRead::NextSegment::FlagBits AlignedRead::NextSegment::compress_flags(const Flags& flags)
{
    FlagBits result {};
    result[0] = flags.is_marked_unmapped;
    result[1] = flags.is_marked_reverse_mapped;
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

AlignedRead::SizeType sequence_size(const AlignedRead& read) noexcept
{
    return static_cast<AlignedRead::SizeType>(read.sequence().size());
}

bool is_empty_sequence(const AlignedRead& read) noexcept
{
    return read.sequence().empty();
}

CigarString splice_cigar(const AlignedRead& read, const GenomicRegion& region)
{
    if (contains(region, read)) return read.cigar_string();
    
    const auto splice_region = overlapped_region(read, region);
    
    return splice(read.cigar_string(), static_cast<CigarOperation::SizeType>(begin_distance(splice_region, read)));
}

AlignedRead::SizeType count_overlapped_bases(const AlignedRead& read, const GenomicRegion& region)
{
    if (contains(region, read)) return sequence_size(read);
    
    // TODO: not quite right as doesn't account for indels
    return static_cast<AlignedRead::SizeType>(std::max(GenomicRegion::DifferenceType {0}, overlap_size(read, region)));
}

bool is_soft_clipped(const AlignedRead& read)
{
    return is_soft_clipped(read.cigar_string());
}

std::pair<AlignedRead::SizeType, AlignedRead::SizeType> get_soft_clipped_sizes(const AlignedRead& read)
{
    return get_soft_clipped_sizes(read.cigar_string());
}

AlignedRead splice(const AlignedRead& read, const GenomicRegion& region)
{
    using std::cbegin;
    
    if (!overlaps(read, region)) {
        throw std::logic_error {"AlignedRead: trying to splice non-overlapped region"};
    }
    
    if (contains(region, read)) return read;
    
    const auto splice_region = overlapped_region(read, region);
    
    const auto reference_offset = static_cast<CigarOperation::SizeType>(begin_distance(splice_region, read));
    
    const auto uncontained_cigar_splice = splice_reference(read.cigar_string(), reference_offset);
    
    auto contained_cigar_splice   = splice_reference(read.cigar_string(), reference_offset,
                                                     region_size(splice_region));
    
    const auto sequence_offset = sequence_size(uncontained_cigar_splice);
    const auto sequence_length = sequence_size(contained_cigar_splice);
    
    AlignedRead::SequenceType sequence_splice(cbegin(read.sequence()) + sequence_offset,
                                              cbegin(read.sequence()) + sequence_offset + sequence_length);
    
    AlignedRead::Qualities qualities_splice(cbegin(read.qualities()) + sequence_offset,
                                            cbegin(read.qualities()) + sequence_offset + sequence_length);
    
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
    return lhs.mapping_quality() == rhs.mapping_quality() &&
           lhs.mapped_region()          == rhs.mapped_region() &&
           lhs.cigar_string()    == rhs.cigar_string() &&
           lhs.qualities()       == rhs.qualities();
}

bool operator<(const AlignedRead& lhs, const AlignedRead& rhs)
{
    if (lhs.mapped_region() == rhs.mapped_region()) {
        if (lhs.mapping_quality() == rhs.mapping_quality()) {
            if (lhs.cigar_string() == rhs.cigar_string()) {
                return lhs.qualities() < rhs.qualities();
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
    if (lhs.is_chimeric() && rhs.is_chimeric()) {
        return lhs.next_segment() == rhs.next_segment();
    }
    return false;
}

bool IsDuplicate::operator()(const AlignedRead &lhs, const AlignedRead &rhs) const
{
    return lhs.mapped_region() == rhs.mapped_region()
        && lhs.cigar_string() == rhs.cigar_string()
        && lhs.flags().is_marked_reverse_mapped == rhs.flags().is_marked_reverse_mapped
        && are_other_segments_duplicates(lhs, rhs);
}

bool operator==(const AlignedRead::NextSegment& lhs, const AlignedRead::NextSegment& rhs)
{
    return lhs.contig_name() == rhs.contig_name()
        && lhs.begin() == rhs.begin()
        && lhs.inferred_template_length() == rhs.inferred_template_length();
}

std::ostream& operator<<(std::ostream& os, const AlignedRead::Qualities& qualities)
{
    std::transform(std::cbegin(qualities), std::cend(qualities),
                   std::ostream_iterator<AlignedRead::QualityType>(os),
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
    if (read.is_chimeric()) {
        os << read.next_segment().contig_name() << ' ';
        os << read.next_segment().begin() << ' ';
        os << read.next_segment().inferred_template_length();
    }
    return os;
}
