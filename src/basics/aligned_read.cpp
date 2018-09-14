// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "aligned_read.hpp"

#include <ostream>

#include <boost/functional/hash.hpp>

#include "utils/sequence_utils.hpp"

namespace octopus {

// AlignedRead::Segment public

const GenomicRegion::ContigName& AlignedRead::Segment::contig_name() const
{
    return contig_name_;
}

GenomicRegion::Position AlignedRead::Segment::begin() const noexcept
{
    return begin_;
}

GenomicRegion::Size AlignedRead::Segment::inferred_template_length() const noexcept
{
    return inferred_template_length_;
}

bool AlignedRead::Segment::is_marked_unmapped() const
{
    return flags_[0];
}

bool AlignedRead::Segment::is_marked_reverse_mapped() const
{
    return flags_[1];
}

// AlignedRead public

const std::string& AlignedRead::name() const noexcept
{
    return name_;
}

const std::string& AlignedRead::read_group() const noexcept
{
    return read_group_;
}

const GenomicRegion& AlignedRead::mapped_region() const noexcept
{
    return region_;
}

const AlignedRead::NucleotideSequence& AlignedRead::sequence() const noexcept
{
    return sequence_;
}

AlignedRead::NucleotideSequence& AlignedRead::sequence() noexcept
{
    return sequence_;
}

const AlignedRead::BaseQualityVector& AlignedRead::base_qualities() const noexcept
{
    return base_qualities_;
}

AlignedRead::BaseQualityVector& AlignedRead::base_qualities() noexcept
{
    return base_qualities_;
}

AlignedRead::MappingQuality AlignedRead::mapping_quality() const noexcept
{
    return mapping_quality_;
}

const CigarString& AlignedRead::cigar() const noexcept
{
    return cigar_;
}

AlignedRead::Direction AlignedRead::direction() const noexcept
{
    if (is_marked_reverse_mapped()) {
        return Direction::reverse;
    } else {
        return Direction::forward;
    }
}

bool AlignedRead::has_other_segment() const noexcept
{
    return static_cast<bool>(next_segment_);
}

const AlignedRead::Segment& AlignedRead::next_segment() const
{
    if (has_other_segment()) {
        return *next_segment_;
    } else {
        throw std::runtime_error {"AlignedRead: read does not have a next segment"};
    }
}

AlignedRead::Flags AlignedRead::flags() const noexcept
{
    return decompress(flags_);
}

void AlignedRead::realign(GenomicRegion new_region, CigarString new_cigar) noexcept
{
    assert(sequence_size(new_cigar) == sequence_.size());
    assert(reference_size(new_cigar) == size(new_region));
    region_ = std::move(new_region);
    cigar_ = std::move(new_cigar);
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

// private methods

AlignedRead::FlagBits AlignedRead::compress(const Flags& flags) const noexcept
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
    result[8] = flags.first_template_segment;
    result[9] = flags.last_template_segment;
    return result;
}

AlignedRead::Flags AlignedRead::decompress(const FlagBits& flags) const noexcept
{
    return {flags[0], flags[1], flags[2], flags[3], flags[4], flags[5], flags[6], flags[7], flags[8], flags[9]};
}

AlignedRead::Segment::FlagBits AlignedRead::Segment::compress(const Flags& flags)
{
    FlagBits result {};
    result[0] = flags.unmapped;
    result[1] = flags.reverse_mapped;
    return result;
}

std::size_t ReadHash::operator()(const octopus::AlignedRead &read) const
{
    std::size_t result {};
    using boost::hash_combine;
    hash_combine(result, std::hash<GenomicRegion>()(read.mapped_region()));
    hash_combine(result, std::hash<CigarString>()(read.cigar()));
    hash_combine(result, boost::hash_range(std::cbegin(read.base_qualities()), std::cend(read.base_qualities())));
    hash_combine(result, read.mapping_quality());
    return result;
}

// Non-member methods

void capitalise_bases(AlignedRead& read) noexcept
{
    utils::capitalise(read.sequence());
}

void cap_qualities(AlignedRead& read, const AlignedRead::BaseQuality max) noexcept
{
    auto& qualities = read.base_qualities();
    std::transform(std::cbegin(qualities), std::cend(qualities), std::begin(qualities),
                   [max] (const auto q) { return std::min(q, max); });
}

void set_front_qualities(AlignedRead& read, std::size_t num_bases, const AlignedRead::BaseQuality value) noexcept
{
    auto& qualities = read.base_qualities();
    std::fill_n(std::begin(qualities), std::min(num_bases, qualities.size()), value);
}

void zero_front_qualities(AlignedRead& read, std::size_t num_bases) noexcept
{
    set_front_qualities(read, num_bases, 0);
}

void set_back_qualities(AlignedRead& read, std::size_t num_bases, const AlignedRead::BaseQuality value) noexcept
{
    auto& qualities = read.base_qualities();
    std::fill_n(std::rbegin(qualities), std::min(num_bases, qualities.size()), value);
}

void zero_back_qualities(AlignedRead& read, std::size_t num_bases) noexcept
{
    set_back_qualities(read, num_bases, 0);
}

bool is_sequence_empty(const AlignedRead& read) noexcept
{
    return read.sequence().empty();
}

AlignedRead::NucleotideSequence::size_type sequence_size(const AlignedRead& read) noexcept
{
    return read.sequence().size();
}

AlignedRead::NucleotideSequence::size_type sequence_size(const AlignedRead& read, const GenomicRegion& region)
{
    if (contig_name(region) != contig_name(read)) return 0;
    if (contains(region, read)) return sequence_size(read);
    const auto copy_region = *overlapped_region(read, region);
    const auto reference_offset = static_cast<CigarOperation::Size>(begin_distance(read, copy_region));
    const auto contained_cigar_copy = copy_reference(read.cigar(), reference_offset, region_size(copy_region));
    return sequence_size(contained_cigar_copy);
}

bool is_forward_strand(const AlignedRead& read) noexcept
{
    return read.direction() == AlignedRead::Direction::forward;
}

bool is_reverse_strand(const AlignedRead& read) noexcept
{
    return !is_forward_strand(read);
}

bool is_soft_clipped(const AlignedRead& read) noexcept
{
    return is_soft_clipped(read.cigar());
}

bool is_front_soft_clipped(const AlignedRead& read) noexcept
{
    return is_front_soft_clipped(read.cigar());
}

bool is_back_soft_clipped(const AlignedRead& read) noexcept
{
    return is_back_soft_clipped(read.cigar());
}

std::pair<CigarOperation::Size, CigarOperation::Size> get_soft_clipped_sizes(const AlignedRead& read) noexcept
{
    return get_soft_clipped_sizes(read.cigar());
}

CigarOperation::Size total_clip_size(const AlignedRead& read) noexcept
{
    const auto p = get_soft_clipped_sizes(read);
    return p.first + p.second;
}

GenomicRegion clipped_mapped_region(const AlignedRead& read)
{
    const auto p = get_soft_clipped_sizes(read);
    using D = GenomicRegion::Distance;
    return expand(mapped_region(read), -static_cast<D>(std::min(p.first, mapped_begin(read))), -static_cast<D>(p.second));
}

bool has_indel(const AlignedRead& read) noexcept
{
    return has_indel(read.cigar());
}

int sum_indel_sizes(const AlignedRead& read) noexcept
{
    return sum_indel_sizes(read.cigar());
}

int max_indel_size(const AlignedRead& read) noexcept
{
    return max_indel_size(read.cigar());
}

CigarString copy_cigar(const AlignedRead& read, const GenomicRegion& region)
{
    if (contains(region, read)) return read.cigar();
    const auto copy_region = *overlapped_region(read, region);
    const auto offset = static_cast<CigarOperation::Size>(begin_distance(read, copy_region));
    return copy(read.cigar(), offset, size(region));
}

AlignedRead copy(const AlignedRead& read, const GenomicRegion& region)
{
    using std::cbegin; using std::next;
    if (!overlaps(read, region)) {
        throw std::logic_error {"AlignedRead: trying to copy non-overlapping region"};
    }
    if (contains(region, read)) return read;
    const auto copy_region = *overlapped_region(read, region);
    const auto reference_offset = static_cast<CigarOperation::Size>(begin_distance(read, copy_region));
    auto uncontained_cigar_copy = copy_reference(read.cigar(), 0, reference_offset);
    auto contained_cigar_copy = copy_reference(read.cigar(), reference_offset, region_size(copy_region));
    if (!uncontained_cigar_copy.empty() && !contained_cigar_copy.empty()
        && uncontained_cigar_copy.back() == contained_cigar_copy.front()
        && is_insertion(uncontained_cigar_copy.back())) {
        uncontained_cigar_copy.pop_back();
    }
    const auto copy_offset = sequence_size(uncontained_cigar_copy);
    const auto copy_length = sequence_size(contained_cigar_copy);
    assert(copy_offset + copy_length <= sequence_size(read));
    const auto subsequence_begin_itr = next(cbegin(read.sequence()), copy_offset);
    const auto subsequence_end_itr   = next(subsequence_begin_itr, copy_length);
    AlignedRead::NucleotideSequence sub_sequence {subsequence_begin_itr, subsequence_end_itr};
    const auto subqualities_begin_itr = next(cbegin(read.base_qualities()), copy_offset);
    const auto subqualities_end_itr   = next(subqualities_begin_itr, copy_length);
    AlignedRead::BaseQualityVector sub_qualities {subqualities_begin_itr, subqualities_end_itr};
    return AlignedRead {read.name(), copy_region, std::move(sub_sequence), std::move(sub_qualities),
                        std::move(contained_cigar_copy), read.mapping_quality(), read.flags(), read.read_group()};
}

template <typename T>
T copy_helper(const T& sequence, const CigarString& cigar, const GenomicRegion& sequence_region, const GenomicRegion& request_region)
{
    if (!overlaps(sequence_region, request_region)) {};
    if (contains(request_region, sequence_region)) return sequence;
    const auto copy_region = *overlapped_region(sequence_region, request_region);
    const auto reference_offset = static_cast<CigarOperation::Size>(begin_distance(sequence_region, copy_region));
    auto uncontained_cigar_copy = copy_reference(cigar, 0, reference_offset);
    auto contained_cigar_copy = copy_reference(cigar, reference_offset, region_size(copy_region));
    if (!uncontained_cigar_copy.empty() && !contained_cigar_copy.empty()
        && uncontained_cigar_copy.back() == contained_cigar_copy.front()
        && is_insertion(uncontained_cigar_copy.back())) {
        uncontained_cigar_copy.pop_back();
    }
    const auto copy_offset = sequence_size(uncontained_cigar_copy);
    const auto copy_length = sequence_size(contained_cigar_copy);
    assert(copy_offset + copy_length <= sequence.size());
    using std::cbegin; using std::next;
    const auto subsequence_begin_itr = next(cbegin(sequence), copy_offset);
    const auto subsequence_end_itr   = next(subsequence_begin_itr, copy_length);
    return {subsequence_begin_itr, subsequence_end_itr};
}

AlignedRead::NucleotideSequence copy_sequence(const AlignedRead& read, const GenomicRegion& region)
{
    return copy_helper(read.sequence(), read.cigar(), read.mapped_region(), region);
}

AlignedRead::BaseQualityVector copy_base_qualities(const AlignedRead& read, const GenomicRegion& region)
{
    return copy_helper(read.base_qualities(), read.cigar(), read.mapped_region(), region);
}

bool operator==(const AlignedRead& lhs, const AlignedRead& rhs) noexcept
{
    return lhs.mapping_quality() == rhs.mapping_quality()
        && lhs.mapped_region()   == rhs.mapped_region()
        && lhs.cigar()           == rhs.cigar()
        && lhs.sequence()        == rhs.sequence()
        && lhs.base_qualities()  == rhs.base_qualities();
}

bool operator<(const AlignedRead& lhs, const AlignedRead& rhs) noexcept
{
    if (lhs.mapped_region() == rhs.mapped_region()) {
        if (lhs.mapping_quality() == rhs.mapping_quality()) {
            if (lhs.cigar() == rhs.cigar()) {
                if (lhs.sequence() == rhs.sequence()) {
                    return lhs.base_qualities() < rhs.base_qualities();
                } else {
                    return lhs.sequence() < rhs.sequence();
                }
            } else {
                return lhs.cigar() < rhs.cigar();
            }
        } else {
            return lhs.mapping_quality() < rhs.mapping_quality();
        }
    } else {
        return lhs.mapped_region() < rhs.mapped_region();
    }
}

bool next_segments_are_duplicates(const AlignedRead& lhs, const AlignedRead& rhs) noexcept
{
    if (lhs.has_other_segment()) {
        if (rhs.has_other_segment()) {
            return lhs.next_segment() == rhs.next_segment();
        } else {
            return false;
        }
    } else {
        return !rhs.has_other_segment();
    }
}

bool IsDuplicate::operator()(const AlignedRead& lhs, const AlignedRead& rhs) const noexcept
{
    return lhs.mapped_region() == rhs.mapped_region()
        && lhs.cigar() == rhs.cigar()
        && lhs.is_marked_reverse_mapped() == rhs.is_marked_reverse_mapped()
        && next_segments_are_duplicates(lhs, rhs);
}

bool operator==(const AlignedRead::Segment& lhs, const AlignedRead::Segment& rhs) noexcept
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
    os << read.name() << ' ';
    os << read.mapped_region() << ' ';
    os << read.sequence() << ' ';
    os << read.base_qualities() << ' ';
    os << read.cigar() << ' ';
    os << static_cast<unsigned>(read.mapping_quality()) << ' ';
    if (read.has_other_segment()) {
        os << read.next_segment().contig_name() << ' ';
        os << read.next_segment().begin() << ' ';
        os << read.next_segment().inferred_template_length();
    }
    return os;
}

} // namespace octopus
