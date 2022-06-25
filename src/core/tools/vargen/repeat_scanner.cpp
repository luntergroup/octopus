// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "repeat_scanner.hpp"

#include <iterator>
#include <algorithm>
#include <numeric>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>

#include "config/common.hpp"
#include "basics/aligned_read.hpp"
#include "basics/cigar_string.hpp"
#include "basics/tandem_repeat.hpp"
#include "io/reference/reference_genome.hpp"
#include "concepts/mappable_range.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/repeat_finder.hpp"
#include "utils/free_memory.hpp"
#include "logging/logging.hpp"

namespace octopus { namespace coretools {

std::unique_ptr<VariantGenerator> RepeatScanner::do_clone() const
{
    return std::make_unique<RepeatScanner>(*this);
}

RepeatScanner::RepeatScanner(const ReferenceGenome& reference, Options options)
: reference_ {reference}
, options_ {options}
, snv_buffer_ {}
, candidates_ {}
{}

bool RepeatScanner::do_requires_reads() const noexcept
{
    return true;
}

void RepeatScanner::do_add_read(const SampleName& sample, const AlignedRead& read)
{
    add_read(read, get_sample_index(sample));
}

void RepeatScanner::do_add_reads(const SampleName& sample, ReadVectorIterator first, ReadVectorIterator last)
{
    const auto sample_index = get_sample_index(sample);
    std::for_each(first, last, [&] (const AlignedRead& read) { add_read(read, sample_index); });
}

void RepeatScanner::do_add_reads(const SampleName& sample, ReadFlatSetIterator first, ReadFlatSetIterator last)
{
    const auto sample_index = get_sample_index(sample);
    std::for_each(first, last, [&] (const AlignedRead& read) { add_read(read, sample_index); });
}

std::vector<Variant> RepeatScanner::do_generate(const RegionSet& regions, OptionalThreadPool workers) const
{
    std::sort(std::begin(candidates_), std::end(candidates_));
    std::vector<Variant> result {};
    for (const auto& region : regions) {
        generate(region, result);
    }
    std::sort(std::begin(result), std::end(result));
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

void RepeatScanner::do_clear() noexcept
{
    free_memory(snv_buffer_);
    free_memory(candidates_);
    free_memory(read_coverage_tracker_);
}

std::string RepeatScanner::name() const
{
    return "RepeatScanner";
}

// private methods

unsigned RepeatScanner::get_sample_index(const SampleName& sample)
{
    const auto itr = std::find(std::cbegin(samples_), std::cend(samples_), sample);
    if (itr != std::cend(samples_)) {
        return std::distance(std::cbegin(samples_), itr);
    } else {
        samples_.push_back(sample);
        return samples_.size() - 1;
    }
}

void RepeatScanner::add_read(const AlignedRead& read, const unsigned sample_index)
{
    assert(snv_buffer_.empty());
    auto ref_index = mapped_begin(read);
    std::size_t read_index {0};
    for (const auto& cigar_operation : read.cigar()) {
        const auto op_size = cigar_operation.size();
        using Flag = CigarOperation::Flag;
        switch (cigar_operation.flag()) {
            case Flag::alignmentMatch:
                add_match_range(GenomicRegion {contig_name(read), ref_index, ref_index + op_size}, read, read_index, sample_index);
                read_index += op_size;
                ref_index  += op_size;
                break;
            case Flag::sequenceMatch:
                read_index += op_size;
                ref_index  += op_size;
                break;
            case Flag::substitution:
            {
                if (op_size == 1) {
                    if (read.base_qualities()[read_index] >= options_.min_base_quality) {
                        add_to_buffer(SNV {ref_index, read.sequence()[read_index]}, contig_name(read), sample_index);
                    }
                } else {
                    reset_buffer(contig_name(read), sample_index);
                }
                read_index += op_size;
                ref_index  += op_size;
                break;
            }
            case Flag::insertion:
            {
                reset_buffer(contig_name(read), sample_index);
                read_index += op_size;
                break;
            }
            case Flag::deletion:
            {
                reset_buffer(contig_name(read), sample_index);
                ref_index += op_size;
                break;
            }
            case Flag::softClipped:
            {
                reset_buffer(contig_name(read), sample_index);
                read_index += op_size;
                ref_index  += op_size;
                break;
            }
            case Flag::hardClipped:
            {
                reset_buffer(contig_name(read), sample_index);
                break;
            }
            case Flag::padding:
                reset_buffer(contig_name(read), sample_index);
                ref_index += op_size;
                break;
            case Flag::skipped:
                reset_buffer(contig_name(read), sample_index);
                ref_index += op_size;
                break;
        }
    }
    reset_buffer(contig_name(read), sample_index);
    read_coverage_tracker_[samples_[sample_index]].add(read);
}

void RepeatScanner::add_match_range(const GenomicRegion& region, const AlignedRead& read, std::size_t read_index, const unsigned sample_index) const
{
    const auto ref_segment = reference_.get().fetch_sequence(region);
    for (std::size_t ref_index {0}; ref_index < ref_segment.size(); ++ref_index, ++read_index) {
        const char ref_base {ref_segment[ref_index]}, read_base {read.sequence()[read_index]};
        if (ref_base != read_base) {
            if (ref_base != 'N' && read_base != 'N') {
                if (read.base_qualities()[read_index] >= options_.min_base_quality) {
                    const auto pos = region.begin() + static_cast<GenomicRegion::Position>(ref_index);
                    add_to_buffer(SNV {pos, read_base}, region.contig_name(), sample_index);
                }
            } else {
                reset_buffer(region.contig_name(), sample_index);
            }
        }
    }
}

void RepeatScanner::add_to_buffer(SNV snv, const ContigName& contig, const unsigned sample_index) const
{
    if (!snv_buffer_.empty()) {
        if (snv_buffer_.front().base != snv.base) {
            reset_buffer(contig, sample_index);
        } else if (snv_buffer_.size() > 1) {
            const auto buffered_snv_gap = begin_distance(snv_buffer_[0], snv_buffer_[1]);
            const auto new_snv_gap = begin_distance(snv_buffer_.back(), snv);
            if (buffered_snv_gap != new_snv_gap) {
                auto tail_snv = snv_buffer_.back();
                reset_buffer(contig, sample_index);
                snv_buffer_.push_back(std::move(tail_snv)); // keep tail snv in case it starts a new run
            }
        } else {
            const auto snv_gap = begin_distance(snv_buffer_.back(), snv);
            if (snv_gap > options_.max_period) {
                reset_buffer(contig, sample_index);
            }
        }
    }
    snv_buffer_.push_back(std::move(snv));
}

void RepeatScanner::reset_buffer(const ContigName& contig, const unsigned sample_index) const
{
    if (snv_buffer_.size() >= options_.min_snvs) {
        assert(!snv_buffer_.empty());
        GenomicRegion region {contig, encompassing_region(snv_buffer_)};
        auto ref_sequence = reference_.get().fetch_sequence(region);
        auto alt_sequence = ref_sequence;
        const auto alt_base = snv_buffer_.front().base;
        alt_sequence[0] = alt_base;
        if (snv_buffer_.size() > 1) {
            const auto period = begin_distance(snv_buffer_[0], snv_buffer_[1]);
            for (std::size_t i {1}; i < snv_buffer_.size(); ++i) {
                alt_sequence[i * period] = alt_base;
            }
        }
        candidates_.emplace_back(Variant {std::move(region), std::move(ref_sequence), std::move(alt_sequence)}, sample_index);
    }
    snv_buffer_.clear();
}

namespace {

struct AdjacentRepeatPair : public Mappable<AdjacentRepeatPair>
{
    TandemRepeat lhs, rhs;
    GenomicRegion region;
    const auto& mapped_region() const noexcept { return region; }
    AdjacentRepeatPair(TandemRepeat lhs, TandemRepeat rhs)
    : lhs {std::move(lhs)}, rhs {std::move(rhs)},
    region {encompassing_region(this->lhs, this->rhs)} {}
};

auto generate_tandem_repeats(const ReferenceGenome& reference, const GenomicRegion& region,
                             const unsigned max_period, const unsigned min_tract_length)
{
    auto result = find_exact_tandem_repeats(reference, region, max_period);
    const auto is_too_short = [=] (const auto& repeat) { return region_size(repeat) < min_tract_length; };
    result.erase(std::remove_if(std::begin(result), std::end(result), is_too_short), std::end(result));
    return result;
}

std::deque<AdjacentRepeatPair>
find_adjacent_tandem_repeats(const ReferenceGenome& reference, const GenomicRegion& region,
                             const unsigned max_period, const unsigned min_tract_length)
{
    const auto repeats = generate_tandem_repeats(reference, region, max_period, min_tract_length);
    std::deque<AdjacentRepeatPair> result {};
    if (repeats.size() > 1) {
        for (auto lhs_itr = std::cbegin(repeats); lhs_itr != std::prev(std::cend(repeats)); ++lhs_itr) {
            const auto possible_adjacents = overlap_range(std::next(lhs_itr), std::cend(repeats), expand_rhs(tail_region(*lhs_itr), 1));
            for (const auto& rhs_repeat : possible_adjacents) {
                if (are_adjacent(*lhs_itr, rhs_repeat)
                    || static_cast<unsigned>(overlap_size(*lhs_itr, rhs_repeat)) < std::max(lhs_itr->period(), rhs_repeat.period())) {
                    result.emplace_back(*lhs_itr, rhs_repeat);
                }
            }
        }
    }
    std::sort(std::begin(result), std::end(result));
    return result;
}

unsigned count_whole_repeats(unsigned mnv_length, unsigned repeat_period) noexcept
{
    return std::ceil(static_cast<double>(mnv_length) / repeat_period);
}

template <typename T>
T repeat(T str, const std::size_t n)
{
    if (n == 0) {
        str.clear();
        str.shrink_to_fit();
        return str;
    } else if (n == 1 || str.empty()) {
        return str;
    }
    const auto period = str.size();
    if (period == 1) {
        str.append(n - 1, str.front());
        return str;
    }
    str.reserve(period * n);
    std::size_t m {2};
    for (; m < n; m *= 2) str += str;
    str.append(str.data(), (n - (m / 2)) * period);
    return str;
}

auto abs_diff(std::size_t lhs, std::size_t rhs) noexcept
{
    return lhs < rhs ? rhs - lhs : lhs - rhs;
}

template <typename Sequence1, typename Sequence2>
bool search(const Sequence1& target, const Sequence2& query)
{
    return std::search(std::cbegin(target), std::cend(target),
                       std::cbegin(query), std::cend(query)) 
                       != std::cend(target);
}

bool can_try_splitting(const Variant& candidate, const AdjacentRepeatPair& repeats)
{
    if (!is_snv(candidate)) return true;
    if (repeats.lhs.period() == 1 || repeats.rhs.period() == 1) {
        if (repeats.lhs.period() > 1) {
            return repeats.rhs.motif() == candidate.alt_allele().sequence()
                && search(repeats.lhs.motif(), candidate.ref_allele().sequence())
                && search(repeats.lhs.motif(), candidate.alt_allele().sequence());
        } else if (repeats.rhs.period() > 1) {
            return repeats.lhs.motif() == candidate.alt_allele().sequence()
                && search(repeats.rhs.motif(), candidate.ref_allele().sequence())
                && search(repeats.rhs.motif(), candidate.alt_allele().sequence());
        } else {
            return search(repeats.lhs.motif(), candidate.alt_allele().sequence())
                || search(repeats.rhs.motif(), candidate.alt_allele().sequence());
        }
    } else {
        return false;
    }
}

} // namespace

void RepeatScanner::generate(const GenomicRegion& region, std::vector<Variant>& result) const
{
    auto mnvs = get_candidate_mnvs(region);
    if (!mnvs.empty()) {
        const auto segments = segment_by_overlapped_move(mnvs);
        for (const auto& segment : segments) {
            assert(!segment.empty());
            const auto segment_region = encompassing_region(segment);
            const auto repeat_search_region = expand(segment_region, 100);
            const auto segment_repeat_pairs = find_adjacent_tandem_repeats(reference_, repeat_search_region, options_.max_period, options_.min_tract_length);
            for (const auto& mnv : segment) {
                for (const auto& repeat_pair : overlap_range(segment_repeat_pairs, mnv)) {
                    if (can_try_splitting(mnv, repeat_pair)) {
                        if (are_adjacent(repeat_pair.lhs, mnv) && contains(repeat_pair.rhs, mnv)) {
                            // insertion of lhs repeat, deletion of rhs repeat
                            const auto num_deleted_periods = count_whole_repeats(region_size(mnv), repeat_pair.rhs.period());
                            auto deleted_region = head_region(repeat_pair.rhs, repeat_pair.rhs.period() * num_deleted_periods);
                            auto deleted_sequence = reference_.get().fetch_sequence(deleted_region);
                            const auto num_inserted_periods = count_whole_repeats(region_size(deleted_region), repeat_pair.lhs.period());
                            auto insertion_region = head_region(repeat_pair.lhs);
                            auto inserted_sequence = repeat(repeat_pair.lhs.motif(), num_inserted_periods);
                            if (deleted_sequence.size() != inserted_sequence.size()) {
                                const auto inbalance = abs_diff(inserted_sequence.size(), deleted_sequence.size());
                                if (repeat_pair.rhs.period() < repeat_pair.lhs.period()) {
                                    if (inbalance % repeat_pair.rhs.period() == 0) {
                                        auto balanced_deletion_region = expand_rhs(deleted_region, inbalance);
                                        auto balanced_deleted_sequence = reference_.get().fetch_sequence(balanced_deletion_region);;
                                        result.emplace_back(std::move(balanced_deletion_region), balanced_deleted_sequence, "");
                                    }
                                } else {
                                    if (inbalance % repeat_pair.lhs.period() == 0) {
                                        auto balanced_insertion_sequence = repeat(repeat_pair.lhs.motif(), num_inserted_periods + inbalance / repeat_pair.lhs.period());
                                        result.emplace_back(insertion_region, "", std::move(balanced_insertion_sequence));
                                    }
                                }
                            }
                            result.emplace_back(std::move(deleted_region), std::move(deleted_sequence), "");
                            result.emplace_back(std::move(insertion_region), "", std::move(inserted_sequence));
                        } else if (are_adjacent(mnv, repeat_pair.rhs) && contains(repeat_pair.lhs, mnv)) {
                            // insertion of rhs repeat, deletion of lhs repeat
                            const auto num_deleted_periods = count_whole_repeats(region_size(mnv), repeat_pair.lhs.period());
                            auto deleted_region = head_region(repeat_pair.lhs, repeat_pair.lhs.period() * num_deleted_periods);
                            auto deleted_sequence = reference_.get().fetch_sequence(deleted_region);
                            const auto num_inserted_periods = count_whole_repeats(region_size(deleted_region), repeat_pair.rhs.period());
                            auto insertion_region = head_region(repeat_pair.rhs);
                            auto inserted_sequence = repeat(repeat_pair.rhs.motif(), num_inserted_periods);
                            if (deleted_sequence.size() != inserted_sequence.size()) {
                                const auto inbalance = abs_diff(inserted_sequence.size(), deleted_sequence.size());
                                if (repeat_pair.rhs.period() < repeat_pair.lhs.period()) {
                                    if (inbalance % repeat_pair.rhs.period() == 0) {
                                        auto balanced_deletion_region = expand_rhs(deleted_region, inbalance);
                                        auto balanced_deleted_sequence = reference_.get().fetch_sequence(balanced_deletion_region);;
                                        result.emplace_back(std::move(balanced_deletion_region), balanced_deleted_sequence, "");
                                    }
                                } else {
                                    if (inbalance % repeat_pair.lhs.period() == 0) {
                                        auto balanced_insertion_sequence = repeat(repeat_pair.lhs.motif(), num_inserted_periods + inbalance / repeat_pair.lhs.period());
                                        result.emplace_back(insertion_region, "", std::move(balanced_insertion_sequence));
                                    }
                                }
                            }
                            result.emplace_back(std::move(deleted_region), std::move(deleted_sequence), "");
                            result.emplace_back(std::move(insertion_region), "", std::move(inserted_sequence));
                        }
                    }
                }
            }
        }
    }
}

std::vector<Variant> RepeatScanner::get_candidate_mnvs(const GenomicRegion& region) const
{
    std::vector<Variant> result {};
    assert(std::is_sorted(std::cbegin(candidates_), std::cend(candidates_)));
    auto viable_candidates = overlap_range(candidates_, region);
    if (empty(viable_candidates)) return result;
    result.reserve(size(viable_candidates, BidirectionallySortedTag {})); // maximum possible
    const auto last_viable_candidate_itr = std::cend(viable_candidates);
    while (!viable_candidates.empty()) {
        const Candidate& candidate {viable_candidates.front()};
        const auto next_candidate_itr = std::find_if_not(std::next(std::cbegin(viable_candidates)), last_viable_candidate_itr,
                                                         [&candidate] (const auto& v) { return v == candidate; });
        const auto num_observations = static_cast<unsigned>(std::distance(std::cbegin(viable_candidates), next_candidate_itr));
        if (num_observations >= options_.min_observations) {
            std::vector<unsigned> sample_observations(samples_.size());
            std::for_each(std::cbegin(viable_candidates), next_candidate_itr, [&] (const auto& v) { ++sample_observations[v.sample_index]; });
            std::vector<unsigned> sample_depths(samples_.size());
            const auto get_depth = [&] (const SampleName& sample) -> unsigned { return read_coverage_tracker_.count(sample) == 1 ? read_coverage_tracker_.at(sample).mean(candidate.mapped_region()) : 0; };
            std::transform(std::cbegin(samples_), std::cend(samples_), std::begin(sample_depths), get_depth);
            std::vector<double> sample_vafs(samples_.size());
            const auto get_vaf = [] (auto count, auto depth) { return static_cast<double>(count) / depth;};
            std::transform(std::cbegin(sample_observations), std::cend(sample_observations), std::cbegin(sample_depths), std::begin(sample_vafs), get_vaf);
            const auto sufficient_sample_observations = [this] (auto obs) { return obs >= options_.min_sample_observations; };
            const auto sufficient_sample_vaf = [this] (auto vaf) { return !options_.min_vaf || vaf >= *options_.min_vaf; };
            if (std::any_of(std::cbegin(sample_observations), std::cend(sample_observations), sufficient_sample_observations)
                && std::any_of(std::cbegin(sample_vafs), std::cend(sample_vafs), sufficient_sample_vaf)) {
                result.push_back(candidate.variant);
            }
        }
        viable_candidates.advance_begin(num_observations);
    }
    return result;
}

} // coretools
} // namespace octopus
