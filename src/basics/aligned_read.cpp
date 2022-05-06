// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "aligned_read.hpp"

#include <ostream>
#include <limits>

#include <boost/functional/hash.hpp>

#include "utils/sequence_utils.hpp"
#include "utils/string_utils.hpp"

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

AlignedRead::Segment::Flags AlignedRead::Segment::flags() const noexcept
{
    return {flags_[0], flags_[1]};
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
    constexpr static Tag readgroup_tag {'R','G'};
    return *this->annotation(readgroup_tag); // RG is required so will always be present
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

CigarString& AlignedRead::cigar() noexcept
{
    return cigar_;
}

AlignedRead::Direction AlignedRead::direction() const noexcept
{
    return is_marked_reverse_mapped() ? Direction::reverse : Direction::forward;
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

void AlignedRead::add_annotation(Tag name, Annotation value)
{
    annoations_.emplace(std::move(name), std::move(value));
}

boost::optional<const AlignedRead::Annotation&> AlignedRead::annotation(const Tag& name) const noexcept
{
    auto itr = annoations_.find(name);
    if (itr != std::cend(annoations_)) {
        return itr->second;
    } else {
        return boost::none;
    }
}

std::vector<AlignedRead::Tag> AlignedRead::tags() const
{
    std::vector<AlignedRead::Tag> result(annoations_.size());
    std::transform(std::cbegin(annoations_), std::cend(annoations_), std::begin(result), 
                   [] (const auto& p) { return p.first; });
    return result;
}

const AlignedRead::NucleotideSequence& AlignedRead::barcode() const noexcept
{
    const static NucleotideSequence empty_barcode {};
    constexpr static Tag barcode_tag {'B','X'};
    auto barcode = this->annotation(barcode_tag);
    return barcode ? *barcode : empty_barcode;
}

void AlignedRead::set_barcode(NucleotideSequence barcode)
{
    constexpr static Tag barcode_tag {'B','X'};
    add_annotation(barcode_tag, std::move(barcode));
}

std::vector<AlignedRead::SupplementaryAlignment> 
AlignedRead::supplementary_alignments() const
{
    constexpr static Tag supplementary_alignment_tag {'S','A'};
    std::vector<AlignedRead::SupplementaryAlignment> result {};
    const auto supplementary_alignments = this->annotation(supplementary_alignment_tag);
    if (supplementary_alignments) {
        for (auto tuple : utils::split(*supplementary_alignments, ';')) {
            auto fields = utils::split(tuple, ',');
            auto cigar = parse_cigar(fields[3]);
            auto begin_pos = boost::lexical_cast<GenomicRegion::Position>(fields[1]);
            if (begin_pos > 0) --begin_pos;
            GenomicRegion region {std::move(fields[0]), begin_pos, begin_pos + reference_size(cigar)};
            auto strand = fields[2] == "+" ? AlignedRead::Direction::forward : AlignedRead::Direction::reverse;
            auto mapping_quality = boost::numeric_cast<AlignedRead::MappingQuality>(boost::lexical_cast<int>(fields[4]));
            result.emplace_back(std::move(region), std::move(cigar), strand, mapping_quality);
        }
    }
    return result;
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

bool AlignedRead::is_marked_next_segment_unmapped() const noexcept
{
    return has_other_segment() && next_segment().is_marked_unmapped();
}

bool AlignedRead::is_marked_reverse_mapped() const noexcept
{
    return flags_[3];
}

bool AlignedRead::is_marked_next_segment_reverse_mapped() const noexcept
{
    return has_other_segment() && next_segment().is_marked_reverse_mapped();
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

bool AlignedRead::is_marked_first_template_segment() const noexcept
{
    return flags_[8];
}

bool AlignedRead::is_marked_last_template_segment() const noexcept
{
    return flags_[9];
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

AlignedRead::SupplementaryAlignment::SupplementaryAlignment(GenomicRegion region, CigarString cigar, AlignedRead::Direction strand, AlignedRead::MappingQuality mapping_quality)
: region_ {std::move(region)}
, cigar_ {std::move(cigar)}
, strand_ {strand}
, mapping_quality_ {mapping_quality}
{}

const GenomicRegion& AlignedRead::SupplementaryAlignment::mapped_region() const noexcept
{
    return region_;
}
const CigarString& AlignedRead::SupplementaryAlignment::cigar() const noexcept
{
    return cigar_;
}
AlignedRead::Direction AlignedRead::SupplementaryAlignment::strand() const noexcept
{
    return strand_;
}
AlignedRead::MappingQuality AlignedRead::SupplementaryAlignment::mapping_quality() const noexcept
{
    return mapping_quality_;
}

// Non-member methods

GenomicRegion::Position five_prime_mapping_position(const AlignedRead& read) noexcept
{
	return is_forward_strand(read) ? mapped_begin(read) : mapped_end(read);
}

GenomicRegion::Position three_prime_mapping_position(const AlignedRead& read) noexcept
{
	return is_reverse_strand(read) ? mapped_begin(read) : mapped_end(read);
}

void capitalise_bases(AlignedRead& read) noexcept
{
    utils::capitalise(read.sequence());
}

unsigned sum_base_qualities(const AlignedRead& read) noexcept
{
	return std::accumulate(std::cbegin(read.base_qualities()), std::cend(read.base_qualities()), 0u);
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
    return read.direction() == AlignedRead::Direction::reverse;
}

bool are_same_strand(const AlignedRead& lhs, const AlignedRead& rhs) noexcept
{
	return lhs.direction() == rhs.direction();
}

bool is_primary_alignment(const AlignedRead& read) noexcept
{
    return !(read.is_marked_secondary_alignment() || read.is_marked_supplementary_alignment());
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
    std::vector<std::pair<AlignedRead::Tag, AlignedRead::Annotation>> annotations {};
    for (const auto& tag : read.tags()) {
        annotations.emplace_back(tag, *read.annotation(tag));
    }
    return AlignedRead {read.name(), copy_region, std::move(sub_sequence), std::move(sub_qualities),
                        std::move(contained_cigar_copy), read.mapping_quality(), read.flags(), read.read_group(),
                        std::move(annotations)};
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

bool is_unlocalized(const AlignedRead::SupplementaryAlignment& alignment) noexcept
{
    return utils::is_suffix("_random", contig_name(alignment));
}

bool is_unplaced(const AlignedRead::SupplementaryAlignment& alignment) noexcept
{
    return utils::is_prefix("chrU_", contig_name(alignment));
}

bool is_decoy(const AlignedRead::SupplementaryAlignment& alignment) noexcept
{
    return utils::is_suffix("_decoy", contig_name(alignment));
}

namespace {

template <typename Range>
Range copy_subsequence(const Range& sequence, std::size_t start, std::size_t len)
{
    auto start_itr = std::next(std::cbegin(sequence), start);
    return {start_itr, std::next(start_itr, std::min(len, sequence.size() - start))};
}

} // namespace

std::vector<AlignedRead> split(const AlignedRead& read, const GenomicRegion::Size chunk_length, const bool append_index_to_read_name)
{
    const auto read_length = sequence_size(read);
    if (chunk_length >= read_length) return {read};
    std::vector<AlignedRead> result {};
    const auto num_chunks = (read_length + chunk_length - 1) / chunk_length;
    result.reserve(num_chunks);
    GenomicRegion::Position chunk_reference_start {mapped_begin(read)};
    for (unsigned idx {0}; idx < num_chunks; ++idx) {
        const auto sequence_offset = idx * chunk_length;
        auto chunk_cigar = copy_sequence(read.cigar(), sequence_offset, chunk_length);
        if (chunk_cigar.empty()) continue;
        // Edge deletions appear at tail and head of adjacent chunks
        if (is_deletion(chunk_cigar.back())) chunk_cigar.pop_back();
        const auto chunk_reference_length = reference_size(chunk_cigar);
        //if (chunk_reference_length == 0) continue;
        auto chunk_sequence = copy_subsequence(read.sequence(), sequence_offset, chunk_length);
        auto chunk_base_quals = copy_subsequence(read.base_qualities(), sequence_offset, chunk_length);
        auto chunk_name = read.name();
        if (append_index_to_read_name) chunk_name += "_" + std::to_string(idx);
        GenomicRegion chunk_region {contig_name(read), chunk_reference_start, chunk_reference_start + chunk_reference_length};
        std::vector<std::pair<AlignedRead::Tag, AlignedRead::Annotation>> annotations {};
        for (auto tag : read.tags()) {
            annotations.emplace_back(tag, *read.annotation(tag));
        }
        if (read.has_other_segment()) {
            result.emplace_back(std::move(chunk_name), std::move(chunk_region), std::move(chunk_sequence),
                                std::move(chunk_base_quals), std::move(chunk_cigar),
                                read.mapping_quality(), read.flags(), read.read_group(),
                                read.next_segment().contig_name(), read.next_segment().begin(),
                                read.next_segment().inferred_template_length(),
                                read.next_segment().flags(), std::move(annotations));
        } else {
            result.emplace_back(std::move(chunk_name), std::move(chunk_region), std::move(chunk_sequence),
                                std::move(chunk_base_quals), std::move(chunk_cigar),
                                read.mapping_quality(), read.flags(), read.read_group(), std::move(annotations));
        }
        chunk_reference_start += chunk_reference_length;
    }
    return result;
}

namespace {

auto calculate_dynamic_bytes(const AlignedRead::Segment& segment)
{
    return segment.contig_name().size() * sizeof(char);
}

auto calculate_dynamic_bytes(const AlignedRead::SupplementaryAlignment& alignment)
{
    return contig_name(alignment).size() * sizeof(char) + alignment.cigar().size() * sizeof(CigarOperation);
}

auto calculate_dynamic_bytes(const std::vector<AlignedRead::SupplementaryAlignment>& alignments)
{
    const static auto add_bytes = [] (auto total, const auto& alignment) { return total + calculate_dynamic_bytes(alignment); };
    return std::accumulate(std::cbegin(alignments), std::cend(alignments), std::size_t {0}, add_bytes);
}

auto calculate_dynamic_bytes(const AlignedRead& read) noexcept
{
    return read.name().size() * sizeof(char)
           + read.read_group().size() * sizeof(char)
           + sequence_size(read) * sizeof(char)
           + sequence_size(read) * sizeof(AlignedRead::BaseQuality)
           + read.cigar().size() * sizeof(CigarOperation)
           + contig_name(read).size() * sizeof(char)
           + read.barcode().size() * sizeof(char)
           + (read.has_other_segment() ? sizeof(AlignedRead::Segment) + calculate_dynamic_bytes(read.next_segment()) : 0)
           + calculate_dynamic_bytes(read.supplementary_alignments());
}

} // namespace

MemoryFootprint footprint(const AlignedRead& read) noexcept
{
    return sizeof(AlignedRead) + calculate_dynamic_bytes(read);
}

bool operator==(const AlignedRead::Segment& lhs, const AlignedRead::Segment& rhs) noexcept
{
    return lhs.contig_name() == rhs.contig_name()
           && lhs.begin() == rhs.begin()
           && lhs.flags_ == rhs.flags_
           && lhs.inferred_template_length() == rhs.inferred_template_length();
}

bool other_segments_equal(const AlignedRead& lhs, const AlignedRead& rhs) noexcept
{
    if (lhs.has_other_segment()) {
        return rhs.has_other_segment() && lhs.next_segment() == rhs.next_segment();
    } else {
        return !rhs.has_other_segment();
    }
}

bool operator==(const AlignedRead& lhs, const AlignedRead& rhs) noexcept
{
    return lhs.mapping_quality() == rhs.mapping_quality()
        && lhs.flags_ == rhs.flags_
        && lhs.mapped_region()   == rhs.mapped_region()
        && lhs.cigar()           == rhs.cigar()
        && lhs.sequence()        == rhs.sequence()
        && lhs.base_qualities()  == rhs.base_qualities()
        && lhs.read_group()      == rhs.read_group()
        && lhs.name()            == rhs.name()
        && other_segments_equal(lhs, rhs);
}

bool operator==(const AlignedRead::Flags& lhs, const AlignedRead::Flags& rhs) noexcept
{
    return lhs.multiple_segment_template    == rhs.multiple_segment_template
        && lhs.all_segments_in_read_aligned == rhs.all_segments_in_read_aligned
        && lhs.unmapped                     == rhs.unmapped
        && lhs.reverse_mapped               == rhs.reverse_mapped
        && lhs.secondary_alignment          == rhs.secondary_alignment
        && lhs.qc_fail                      == rhs.qc_fail
        && lhs.duplicate                    == rhs.duplicate
        && lhs.supplementary_alignment      == rhs.supplementary_alignment
        && lhs.first_template_segment       == rhs.first_template_segment
        && lhs.last_template_segment        == rhs.last_template_segment;
        
}

bool operator<(const AlignedRead& lhs, const AlignedRead& rhs) noexcept
{
    if (lhs.mapped_region() == rhs.mapped_region()) {
        if (lhs.direction() == rhs.direction()) {
            if (lhs.mapping_quality() == rhs.mapping_quality()) {
                if (lhs.cigar() == rhs.cigar()) {
                    if (lhs.sequence() == rhs.sequence()) {
                        if (lhs.read_group() == rhs.read_group()) {
                            if (lhs.name() == rhs.name()) {
                                return lhs.base_qualities() < rhs.base_qualities();
                            } else {
                                return lhs.name() < rhs.name();
                            }
                        } else {
                            return lhs.read_group() < rhs.read_group();
                        }
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
            return lhs.direction() == AlignedRead::Direction::forward; // put forward strand reads first
        }
    } else {
        return lhs.mapped_region() < rhs.mapped_region();
    }
}

std::ostream& operator<<(std::ostream& os, const AlignedRead::BaseQualityVector& qualities)
{
    std::transform(std::cbegin(qualities), std::cend(qualities),
                   std::ostream_iterator<AlignedRead::BaseQuality>(os),
                   [] (const auto q) { return static_cast<unsigned>(q + 33); }
                   );
    return os;
}

namespace {

/*! @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
static constexpr std::uint16_t BAM_FPAIRED =       1;
/*! @abstract the read is mapped in a proper pair */
static constexpr std::uint16_t BAM_FPROPER_PAIR =  2;
/*! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
static constexpr std::uint16_t BAM_FUNMAP =        4;
/*! @abstract the mate is unmapped */
static constexpr std::uint16_t BAM_FMUNMAP =       8;
/*! @abstract the read is mapped to the reverse strand */
static constexpr std::uint16_t BAM_FREVERSE =     16;
/*! @abstract the mate is mapped to the reverse strand */
static constexpr std::uint16_t BAM_FMREVERSE =    32;
/*! @abstract this is read1 */
static constexpr std::uint16_t BAM_FREAD1 =       64;
/*! @abstract this is read2 */
static constexpr std::uint16_t BAM_FREAD2 =      128;
/*! @abstract not primary alignment */
static constexpr std::uint16_t BAM_FSECONDARY =  256;
/*! @abstract QC failure */
static constexpr std::uint16_t BAM_FQCFAIL =     512;
/*! @abstract optical or PCR duplicate */
static constexpr std::uint16_t BAM_FDUP =       1024;
/*! @abstract supplementary alignment */
static constexpr std::uint16_t BAM_FSUPPLEMENTARY = 2048;

void set_flag(bool set, std::uint16_t mask, std::uint16_t& result) noexcept
{
    constexpr std::uint16_t zeros {0}, ones {std::numeric_limits<std::uint16_t>::max()};
    result |= (set ? ones : zeros) & mask;
}

auto compute_flag_bits(const AlignedRead& read) noexcept
{
    std::uint16_t result {};
    set_flag(read.is_marked_multiple_segment_template(),    BAM_FPAIRED, result);
    set_flag(read.is_marked_all_segments_in_read_aligned(), BAM_FPROPER_PAIR, result);
    set_flag(read.is_marked_unmapped(),                     BAM_FUNMAP, result);
    set_flag(read.is_marked_next_segment_unmapped(),        BAM_FMUNMAP, result);
    set_flag(read.is_marked_reverse_mapped(),               BAM_FREVERSE, result);
    set_flag(read.is_marked_next_segment_reverse_mapped(),  BAM_FMREVERSE, result);
    set_flag(read.is_marked_secondary_alignment(),          BAM_FSECONDARY, result);
    set_flag(read.is_marked_qc_fail(),                      BAM_FQCFAIL, result);
    set_flag(read.is_marked_duplicate(),                    BAM_FDUP, result);
    set_flag(read.is_marked_supplementary_alignment(),      BAM_FSUPPLEMENTARY, result);
    set_flag(read.is_marked_first_template_segment(),       BAM_FREAD1, result);
    set_flag(read.is_marked_last_template_segment(),        BAM_FREAD2, result);
    return result;
}

} // namespace

std::ostream& operator<<(std::ostream& os, const AlignedRead& read)
{
    os << read.name() << '\t';
    os << compute_flag_bits(read) << '\t';
    os << contig_name(read) << '\t';
    os << mapped_begin(read) << '\t';
    os << static_cast<unsigned>(read.mapping_quality()) << '\t';
    os << read.cigar() << '\t';
    if (read.has_other_segment()) {
        os << read.next_segment().contig_name() << '\t';
        os << read.next_segment().begin() << '\t';
        os << read.next_segment().inferred_template_length() << '\t';
    } else {
        os << "*\t0\t0\t";
    }
    os << read.sequence() << '\t';
    os << read.base_qualities() << '\t';
    return os;
}

} // namespace octopus
