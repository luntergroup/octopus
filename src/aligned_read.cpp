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

const std::string& AlignedRead::NextSegment::get_contig_name() const
{
    return contig_name_;
}

AlignedRead::SizeType AlignedRead::NextSegment::get_begin() const noexcept
{
    return begin_;
}

AlignedRead::SizeType AlignedRead::NextSegment::get_inferred_template_length() const noexcept
{
    return inferred_template_length_;
}

GenomicRegion AlignedRead::NextSegment::get_inferred_region() const
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

const GenomicRegion& AlignedRead::get_region() const noexcept
{
    return region_;
}

const AlignedRead::SequenceType& AlignedRead::get_sequence() const noexcept
{
    return sequence_;
}

const AlignedRead::Qualities& AlignedRead::get_qualities() const noexcept
{
    return qualities_;
}

AlignedRead::QualityType AlignedRead::get_mapping_quality() const noexcept
{
    return mapping_quality_;
}

AlignedRead::SizeType AlignedRead::get_sequence_size() const noexcept
{
    return static_cast<SizeType>(sequence_.size());
}

const CigarString& AlignedRead::get_cigar_string() const noexcept
{
    return cigar_string_;
}

bool AlignedRead::has_mate() const noexcept
{
    return next_segment_ != nullptr;
}

const AlignedRead::NextSegment& AlignedRead::get_next_segment() const
{
    if (is_chimeric()) {
        return *next_segment_;
    } else {
        throw std::runtime_error {"Read does not have a next segment"};
    }
}

AlignedRead::Flags AlignedRead::get_flags() const
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

bool AlignedRead::is_chimeric() const noexcept
{
    return next_segment_ != nullptr;
}

bool AlignedRead::is_marked_all_segments_in_read_aligned() const
{
    return flags_[0];
}

bool AlignedRead::is_marked_multiple_read_template() const
{
    return flags_[1];
}

bool AlignedRead::is_marked_unmapped() const
{
    return flags_[2];
}

bool AlignedRead::is_marked_reverse_mapped() const
{
    return flags_[3];
}

bool AlignedRead::is_marked_secondary_alignment() const
{
    return flags_[4];
}

bool AlignedRead::is_marked_qc_fail() const
{
    return flags_[5];
}

bool AlignedRead::is_marked_duplicate() const
{
    return flags_[6];
}

bool AlignedRead::is_marked_supplementary_alignment() const
{
    return flags_[7];
}

size_t AlignedRead::get_hash() const
{
    if (hash_ == 0) { // 0 is reserved
        hash_ = make_hash(); // lazy evaluation
    }
    return hash_;
}

void AlignedRead::zero_front_qualities(const SizeType num_bases) noexcept
{
    std::for_each(std::begin(qualities_),
                  std::next(std::begin(qualities_), std::min(num_bases, get_sequence_size())),
                  [] (auto& quality) { quality = 0; });
}

void AlignedRead::zero_back_qualities(const SizeType num_bases) noexcept
{
    std::for_each(std::rbegin(qualities_),
                  std::next(std::rbegin(qualities_), std::min(num_bases, get_sequence_size())),
                  [] (auto& quality) { quality = 0; });
}

// TODO: can we also compress qualities and cigar string?
void AlignedRead::compress()
{
    sequence_ = Octopus::compress(sequence_);
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

size_t AlignedRead::make_hash() const
{
    using boost::hash_combine;
    
    size_t result {};
    
    hash_combine(result, std::hash<GenomicRegion>()(region_));
    hash_combine(result, std::hash<CigarString>()(cigar_string_));
    hash_combine(result, boost::hash_range(std::cbegin(qualities_), std::cend(qualities_)));
    hash_combine(result, mapping_quality_);
    
    if (result != 0) return result;
    
    hash_combine(result, std::hash<CigarString>()(cigar_string_));
    
    return (result == 0) ? 1 : result; // 0 is reserved
}

// Non-member methods

CigarString splice_cigar(const AlignedRead& read, const GenomicRegion& region)
{
    if (contains(region, read)) return read.get_cigar_string();
    
    const auto splice_region = get_overlapped(read, region);
    
    return splice(read.get_cigar_string(), get_begin(splice_region) - get_begin(read), size(splice_region));
}

AlignedRead::SizeType count_overlapped_bases(const AlignedRead& read, const GenomicRegion& region)
{
    if (contains(region, read)) return read.get_sequence_size();
    
    // TODO: not quite right as doesn't account for indels
    return static_cast<AlignedRead::SizeType>(std::max(GenomicRegion::DifferenceType {0}, overlap_size(read, region)));
}

bool is_soft_clipped(const AlignedRead& read)
{
    return is_soft_clipped(read.get_cigar_string());
}

std::pair<AlignedRead::SizeType, AlignedRead::SizeType> get_soft_clipped_sizes(const AlignedRead& read)
{
    return get_soft_clipped_sizes(read.get_cigar_string());
}

AlignedRead splice(const AlignedRead& read, const GenomicRegion& region)
{
    using std::cbegin;
    
    if (!overlaps(read, region)) {
        throw std::logic_error {"AlignedRead: trying to splice non-overlapped region"};
    }
    
    if (contains(region, read)) return read;
    
    const auto splice_region = get_overlapped(read, region);
    
    const auto reference_offset = get_begin(splice_region) - get_begin(read);
    
    const auto uncontained_cigar_splice = splice_reference(read.get_cigar_string(), reference_offset);
    
    auto contained_cigar_splice   = splice_reference(read.get_cigar_string(), reference_offset,
                                                     size(splice_region));
    
    const auto sequence_offset = sequence_size(uncontained_cigar_splice);
    const auto sequence_length = sequence_size(contained_cigar_splice);
    
    AlignedRead::SequenceType sequence_splice(cbegin(read.get_sequence()) + sequence_offset,
                                              cbegin(read.get_sequence()) + sequence_offset + sequence_length);
    
    AlignedRead::Qualities qualities_splice(cbegin(read.get_qualities()) + sequence_offset,
                                            cbegin(read.get_qualities()) + sequence_offset + sequence_length);
    
    return AlignedRead {
        splice_region,
        std::move(sequence_splice),
        std::move(qualities_splice),
        std::move(contained_cigar_splice),
        read.get_mapping_quality(),
        read.get_flags()
    };
}

bool operator==(const AlignedRead& lhs, const AlignedRead& rhs)
{
    return lhs.get_mapping_quality() == rhs.get_mapping_quality() &&
           lhs.get_region()          == rhs.get_region() &&
           lhs.get_cigar_string()    == rhs.get_cigar_string() &&
           lhs.get_qualities()       == rhs.get_qualities();
}

bool operator<(const AlignedRead& lhs, const AlignedRead& rhs)
{
    if (lhs.get_region() == rhs.get_region()) {
        if (lhs.get_mapping_quality() == rhs.get_mapping_quality()) {
            if (lhs.get_cigar_string() == rhs.get_cigar_string()) {
                return lhs.get_qualities() < rhs.get_qualities();
            } else {
                return lhs.get_cigar_string() < rhs.get_cigar_string();
            }
        } else {
            return lhs.get_mapping_quality() < rhs.get_mapping_quality();
        }
    } else {
        return lhs.get_region() < rhs.get_region();
    }
}

bool IsDuplicate::operator()(const AlignedRead &lhs, const AlignedRead &rhs) const
{
    return lhs.get_region() == rhs.get_region() && lhs.get_cigar_string() == rhs.get_cigar_string();
}

bool operator==(const AlignedRead::NextSegment& lhs, const AlignedRead::NextSegment& rhs)
{
    return lhs.get_contig_name() == rhs.get_contig_name() && lhs.get_begin() == rhs.get_begin();
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
    os << read.get_region() << '\n';
    os << read.get_sequence() << '\n';
    os << read.get_qualities() << '\n';
    os << read.get_cigar_string() << '\n';
    os << static_cast<unsigned>(read.get_mapping_quality()) << '\n';
    if (read.is_chimeric()) {
        os << read.get_next_segment().get_contig_name() << '\n';
        os << read.get_next_segment().get_begin() << '\n';
        os << read.get_next_segment().get_inferred_template_length();
    } else {
        os << "no other segments";
    }
    return os;
}
