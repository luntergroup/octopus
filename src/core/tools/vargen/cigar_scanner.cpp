// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "cigar_scanner.hpp"

#include <iterator>
#include <algorithm>
#include <numeric>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>

#include "config/common.hpp"
#include "basics/aligned_read.hpp"
#include "basics/cigar_string.hpp"
#include "io/reference/reference_genome.hpp"
#include "concepts/mappable_range.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/append.hpp"
#include "utils/sequence_utils.hpp"
#include "utils/free_memory.hpp"
#include "logging/logging.hpp"

#include "utils/maths.hpp"

namespace octopus { namespace coretools {

std::unique_ptr<VariantGenerator> CigarScanner::do_clone() const
{
    return std::make_unique<CigarScanner>(*this);
}

CigarScanner::CigarScanner(const ReferenceGenome& reference, Options options)
: reference_ {reference}
, options_ {options}
, buffer_ {}
, candidates_ {}
, likely_misaligned_candidates_ {}
, max_seen_candidate_size_ {}
, combined_read_coverage_tracker_ {}
, misaligned_read_coverage_tracker_ {}
, sample_read_coverage_tracker_ {}
, sample_forward_strand_coverage_tracker_ {}
, artificial_read_buffer_ {}
{
    buffer_.reserve(100);
}

bool CigarScanner::do_requires_reads() const noexcept
{
    return true;
}

namespace {

template <typename Sequence, typename P, typename S>
Sequence copy(const Sequence& sequence, const P pos, const S size)
{
    const auto it = std::next(std::cbegin(sequence), pos);
    return Sequence {it, std::next(it, size)};
}

double ln_probability_read_correctly_aligned(const double misalign_penalty, const AlignedRead& read,
                                             const double max_expected_mutation_rate)
{
    const auto k = static_cast<unsigned>(std::floor(misalign_penalty));
    if (k == 0) {
        return 0;
    } else {
        const auto ln_prob_missmapped = -maths::constants::ln10Div10<> * read.mapping_quality();
        const auto ln_prob_mapped = std::log(1.0 - std::exp(ln_prob_missmapped));
        const auto mu = max_expected_mutation_rate * region_size(read);
        auto ln_prob_given_mapped = maths::log_poisson_sf(k, mu);
        return ln_prob_mapped + ln_prob_given_mapped;
    }
}

} // namespace

void CigarScanner::do_add_read(const SampleName& sample, const AlignedRead& read)
{
    add_read(sample, read, sample_read_coverage_tracker_[sample], sample_forward_strand_coverage_tracker_[sample]);
}
void CigarScanner::do_add_template(const SampleName& sample, const AlignedTemplate& reads)
{
    add_template(sample, reads, sample_read_coverage_tracker_[sample], sample_forward_strand_coverage_tracker_[sample]);
}

void CigarScanner::add_read(const SampleName& sample, const AlignedRead& read,
                            CoverageTracker<GenomicRegion>& coverage_tracker,
                            CoverageTracker<GenomicRegion>& forward_strand_coverage_tracker)
{
    using std::cbegin; using std::next; using std::move;
    using Flag = CigarOperation::Flag;
    const auto& read_contig   = contig_name(read);
    const auto& read_sequence = read.sequence();
    auto ref_index = mapped_begin(read);
    std::size_t read_index {0};
    GenomicRegion region;
    double misalignment_penalty {0};
    buffer_.clear();
    for (const auto& cigar_operation : read.cigar()) {
        const auto op_size = cigar_operation.size();
        switch (cigar_operation.flag()) {
            case Flag::alignmentMatch:
                misalignment_penalty += add_snvs_in_match_range(GenomicRegion {read_contig, ref_index, ref_index + op_size},
                                                                read, read_index, sample);
                read_index += op_size;
                ref_index  += op_size;
                break;
            case Flag::sequenceMatch:
                read_index += op_size;
                ref_index  += op_size;
                break;
            case Flag::substitution:
            {
                region = GenomicRegion {read_contig, ref_index, ref_index + op_size};
                auto ref_sequence = reference_.get().fetch_sequence(region);
                if (ref_sequence.size() > 1 && options_.split_mnvs) {
                    for (CigarOperation::Size snv_offset {0}; snv_offset < op_size; ++snv_offset) {
                        add_candidate(GenomicRegion {read_contig, ref_index + snv_offset, ref_index + snv_offset + 1},
                                      ref_sequence[snv_offset],
                                      copy(read_sequence, read_index + snv_offset, 1),
                                      read, read_index + snv_offset, sample);
                    }
                } else {
                    add_candidate(std::move(region),
                                  std::move(ref_sequence),
                                  copy(read_sequence, read_index, op_size),
                                  read, read_index, sample);
                }
                read_index += op_size;
                ref_index  += op_size;
                if (options_.misalignment_parameters) misalignment_penalty += op_size * options_.misalignment_parameters->snv_penalty;
                break;
            }
            case Flag::insertion:
            {
                add_candidate(GenomicRegion {read_contig, ref_index, ref_index},
                              "",
                              copy(read_sequence, read_index, op_size),
                              read, read_index, sample);
                read_index += op_size;
                if (options_.misalignment_parameters)  misalignment_penalty += options_.misalignment_parameters->indel_penalty;
                break;
            }
            case Flag::deletion:
            {
                region = GenomicRegion {read_contig, ref_index, ref_index + op_size};
                add_candidate(move(region),
                              reference_.get().fetch_sequence(region),
                              "",
                              read, read_index, sample);
                ref_index += op_size;
                if (options_.misalignment_parameters)  misalignment_penalty += options_.misalignment_parameters->indel_penalty;
                break;
            }
            case Flag::softClipped:
            {
                read_index += op_size;
                ref_index  += op_size;
                if (options_.misalignment_parameters && op_size > options_.misalignment_parameters->max_unpenalised_clip_size) {
                    misalignment_penalty += options_.misalignment_parameters->clip_penalty;
                }
                break;
            }
            case Flag::hardClipped:
            {
                if (options_.misalignment_parameters && op_size > options_.misalignment_parameters->max_unpenalised_clip_size) {
                    misalignment_penalty += options_.misalignment_parameters->clip_penalty;
                }
                break;
            }
            case Flag::padding:
                ref_index += op_size;
                break;
            case Flag::skipped:
                ref_index += op_size;
                break;
        }
    }
    if (options_.use_clipped_coverage_tracking) {
        const auto clipped_region = clipped_mapped_region(read);
        combined_read_coverage_tracker_.add(clipped_region);
        coverage_tracker.add(clipped_region);
        if (is_forward_strand(read)) forward_strand_coverage_tracker.add(clipped_region);
    } else {
        combined_read_coverage_tracker_.add(read);
        coverage_tracker.add(read);
        if (is_forward_strand(read)) forward_strand_coverage_tracker.add(read);
    }
    if (!is_likely_misaligned(read, misalignment_penalty)) {
        utils::append(std::move(buffer_), candidates_);
    } else {
        utils::append(std::move(buffer_), likely_misaligned_candidates_);
        misaligned_read_coverage_tracker_.add(clipped_mapped_region(read));
    }
}

bool have_split_insertion(const AlignedRead& lhs, const AlignedRead& rhs)
{
    return are_adjacent(lhs, rhs) && is_insertion(lhs.cigar().back()) && is_insertion(rhs.cigar().front()) && rhs.cigar().size() == 1;
}

void CigarScanner::add_template(const SampleName& sample, const AlignedTemplate& reads,
                                CoverageTracker<GenomicRegion>& coverage_tracker,
                                CoverageTracker<GenomicRegion>& forward_strand_coverage_tracker)
{
    for (auto read_itr = std::cbegin(reads); read_itr != std::cend(reads);) {
        const static auto has_tail_insertion = [] (const AlignedRead& read) { return is_insertion(read.cigar().back()); };
        const static auto is_all_insertion = [] (const AlignedRead& read) { return read.cigar().size() == 1 && is_insertion(read.cigar().front()); };
        const static auto has_head_insertion = [] (const AlignedRead& read) { return is_insertion(read.cigar().front()); };
        const auto split_insertion_begin_itr = std::find_if(read_itr, std::cend(reads), has_tail_insertion);
        std::for_each(read_itr, split_insertion_begin_itr, [&] (const auto& read) {
            add_read(sample, read, coverage_tracker, forward_strand_coverage_tracker);
        });
        if (split_insertion_begin_itr != std::cend(reads)) {
            auto split_insertion_end_itr = std::find_if_not(std::next(split_insertion_begin_itr), std::cend(reads), is_all_insertion);   
            AlignedRead extended_read {*split_insertion_begin_itr};
            std::for_each(std::next(split_insertion_begin_itr), split_insertion_end_itr, [&] (const auto& read) {
                assert(read.cigar().size() == 1);
                utils::append(read.sequence(), extended_read.sequence());
                utils::append(read.base_qualities(), extended_read.base_qualities());
                extended_read.cigar().back().set_size(extended_read.cigar().back().size() + read.cigar().front().size());
            });
            if (split_insertion_end_itr != std::cend(reads) && has_head_insertion(*split_insertion_end_itr)) {
                AlignedRead chopped_read {*split_insertion_end_itr};
                ++split_insertion_end_itr;
                const auto chop_length = chopped_read.cigar().front().size();
                extended_read.cigar().back().set_size(extended_read.cigar().back().size() + chop_length);
                chopped_read.cigar().erase(std::cbegin(chopped_read.cigar()));
                extended_read.sequence().insert(std::cend(extended_read.sequence()), std::cbegin(chopped_read.sequence()), std::next(std::cbegin(chopped_read.sequence()), chop_length));
                chopped_read.sequence().erase(std::cbegin(chopped_read.sequence()), std::next(std::cbegin(chopped_read.sequence()), chop_length));
                extended_read.base_qualities().insert(std::cend(extended_read.base_qualities()), std::cbegin(chopped_read.base_qualities()), std::next(std::cbegin(chopped_read.base_qualities()), chop_length));
                chopped_read.base_qualities().erase(std::cbegin(chopped_read.base_qualities()), std::next(std::cbegin(chopped_read.base_qualities()), chop_length));
                artificial_read_buffer_.push_back(std::move(extended_read));
                add_read(sample, artificial_read_buffer_.back(), coverage_tracker, forward_strand_coverage_tracker);
                artificial_read_buffer_.push_back(std::move(chopped_read));
                add_read(sample, artificial_read_buffer_.back(), coverage_tracker, forward_strand_coverage_tracker);
            } else {
                artificial_read_buffer_.push_back(std::move(extended_read));
                add_read(sample, artificial_read_buffer_.back(), coverage_tracker, forward_strand_coverage_tracker);
            }
            read_itr = split_insertion_end_itr;
        } else {
            break;
        }
    }
}

void CigarScanner::do_add_reads(const SampleName& sample, ReadVectorIterator first, ReadVectorIterator last)
{
    auto& coverage_tracker = sample_read_coverage_tracker_[sample];
    auto& forward_strand_coverage_tracker = sample_forward_strand_coverage_tracker_[sample];
    std::for_each(first, last, [&] (const auto& read) { add_read(sample, read, coverage_tracker, forward_strand_coverage_tracker); });
}
void CigarScanner::do_add_reads(const SampleName& sample, ReadFlatSetIterator first, ReadFlatSetIterator last)
{
    auto& coverage_tracker = sample_read_coverage_tracker_[sample];
    auto& forward_strand_coverage_tracker = sample_forward_strand_coverage_tracker_[sample];
    std::for_each(first, last, [&] (const auto& read) { add_read(sample, read, coverage_tracker, forward_strand_coverage_tracker); });
}
void CigarScanner::do_add_reads(const SampleName& sample, TemplateVectorIterator first, TemplateVectorIterator last)
{
    auto& coverage_tracker = sample_read_coverage_tracker_[sample];
    auto& forward_strand_coverage_tracker = sample_forward_strand_coverage_tracker_[sample];
    std::for_each(first, last, [&] (const auto& reads) { add_template(sample, reads, coverage_tracker, forward_strand_coverage_tracker); });
}
void CigarScanner::do_add_reads(const SampleName& sample, TemplateFlatSetIterator first, TemplateFlatSetIterator last)
{
    auto& coverage_tracker = sample_read_coverage_tracker_[sample];
    auto& forward_strand_coverage_tracker = sample_forward_strand_coverage_tracker_[sample];
    std::for_each(first, last, [&] (const auto& reads) { add_template(sample, reads, coverage_tracker, forward_strand_coverage_tracker); });
}

unsigned get_min_depth(const Variant& v, const CoverageTracker<GenomicRegion>& tracker)
{
    if (is_insertion(v)) {
        const auto& region = mapped_region(v);
        if (region.begin() > 0) {
            return tracker.min(expand(region, 1, 1));
        } else {
            return tracker.min(expand_rhs(region, 1));
        }
    } else {
        return tracker.min(mapped_region(v));
    }
}

struct VariantBucket : public Mappable<VariantBucket>
{
    VariantBucket(GenomicRegion region) : region {std::move(region)} {}
    GenomicRegion region;
    std::deque<Variant> variants;
    const GenomicRegion& mapped_region() const noexcept { return region; }
};

auto init_variant_buckets(const std::vector<GenomicRegion>& regions)
{
    std::vector<VariantBucket> result {};
    result.reserve(regions.size());
    std::transform(std::cbegin(regions), std::cend(regions), std::back_inserter(result),
                   [] (const auto& region) { return VariantBucket {region}; });
    return result;
}

boost::optional<VariantBucket&> find_contained(std::vector<VariantBucket>& buckets, const Variant& variant)
{
    if (buckets.empty()) return boost::none;
    const auto itr = std::find_if(std::begin(buckets), std::end(buckets),
                                  [&] (const auto& region) { return contains(region, variant); });
    if (itr != std::end(buckets)) {
        assert(contains(*itr, variant));
        return *itr;
    } else {
        return boost::none;
    }
}

void choose_push_back(Variant candidate, std::vector<Variant>& final_candidates,
                      std::vector<VariantBucket>& repeat_buckets)
{
    auto bucket = find_contained(repeat_buckets, candidate);
    if (bucket) {
        bucket->variants.push_back(std::move(candidate));
    } else {
        final_candidates.push_back(std::move(candidate));
    }
}

std::vector<Variant> CigarScanner::do_generate(const RegionSet& regions, OptionalThreadPool workers) const
{
    std::sort(std::begin(candidates_), std::end(candidates_));
    std::sort(std::begin(likely_misaligned_candidates_), std::end(likely_misaligned_candidates_));
    std::vector<Variant> result {};
    for (const auto& region : regions) {
        generate(region, result);
    }
    return result;
}

void CigarScanner::do_clear() noexcept
{
    free_memory(buffer_);
    free_memory(candidates_);
    free_memory(likely_misaligned_candidates_);
    free_memory(combined_read_coverage_tracker_);
    free_memory(misaligned_read_coverage_tracker_);
    free_memory(sample_read_coverage_tracker_);
    free_memory(sample_forward_strand_coverage_tracker_);
    free_memory(artificial_read_buffer_);
    max_seen_candidate_size_ = 0;
}

std::string CigarScanner::name() const
{
    return "CigarScanner";
}

// private methods

double CigarScanner::add_snvs_in_match_range(const GenomicRegion& region, const AlignedRead& read,
                                             std::size_t read_index, const SampleName& origin)
{
    const NucleotideSequence ref_segment {reference_.get().fetch_sequence(region)};
    double misalignment_penalty {0};
    for (std::size_t ref_index {0}; ref_index < ref_segment.size(); ++ref_index, ++read_index) {
        const char ref_base {ref_segment[ref_index]}, read_base {read.sequence()[read_index]};
        if (ref_base != read_base && ref_base != 'N' && read_base != 'N') {
            const auto begin_pos = region.begin() + static_cast<GenomicRegion::Position>(ref_index);
            add_candidate(GenomicRegion {region.contig_name(), begin_pos, begin_pos + 1},
                          ref_base, read_base, read, read_index, origin);
            if (options_.misalignment_parameters && read.base_qualities()[read_index] >= options_.misalignment_parameters->snv_threshold) {
                misalignment_penalty += options_.misalignment_parameters->snv_penalty;
            }
        }
    }
    return misalignment_penalty;
}

void CigarScanner::generate(const GenomicRegion& region, std::vector<Variant>& result) const
{
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::next;
    assert(std::is_sorted(std::cbegin(candidates_), std::cend(candidates_)));
    auto viable_candidates = overlap_range(candidates_, region, max_seen_candidate_size_);
    if (empty(viable_candidates)) return;
    result.reserve(result.size() + size(viable_candidates, BidirectionallySortedTag {})); // maximum possible
    const auto last_viable_candidate_itr = cend(viable_candidates);
    while (!viable_candidates.empty()) {
        const Candidate& candidate {viable_candidates.front()};
        const auto next_candidate_itr = std::find_if_not(next(cbegin(viable_candidates)), last_viable_candidate_itr,
                                                         [this, &candidate] (const Candidate& c) {
                                                             return options_.match(c.variant, candidate.variant);
                                                         });
        const auto num_matches = std::distance(cbegin(viable_candidates), next_candidate_itr);
        const auto observation = make_observation(cbegin(viable_candidates), next_candidate_itr);
        if (options_.include(observation)) {
            if (num_matches > 1) {
                auto unique_itr = cbegin(viable_candidates);
                while (unique_itr != next_candidate_itr) {
                    result.push_back(unique_itr->variant);
                    unique_itr = std::find_if_not(next(unique_itr), next_candidate_itr,
                                                  [unique_itr] (const Candidate& c) {
                                                      return c.variant == unique_itr->variant;
                                                  });
                }
            } else {
                result.push_back(candidate.variant);
            }
        }
        viable_candidates.advance_begin(num_matches);
    }
    if (debug_log_ && !likely_misaligned_candidates_.empty()) {
        const auto novel_unique_misaligned_variants = get_novel_likely_misaligned_candidates(result);
        if (!novel_unique_misaligned_variants.empty()) {
            stream(*debug_log_) << "DynamicCigarScanner: ignoring "
                                << count_overlapped(novel_unique_misaligned_variants, region)
                                << " unique candidates in " << region;
        }
    }
}

unsigned CigarScanner::sum_base_qualities(const Candidate& candidate) const noexcept
{
    const auto first_base_qual_itr = std::next(std::cbegin(candidate.source.get().base_qualities()), candidate.offset);
    const auto last_base_qual_itr = std::next(first_base_qual_itr, alt_sequence_size(candidate.variant));
    return std::accumulate(first_base_qual_itr, last_base_qual_itr, 0u);
}

bool CigarScanner::is_likely_misaligned(const AlignedRead& read, const double penalty) const
{
    if (options_.misalignment_parameters) {
        auto mu = options_.misalignment_parameters->max_expected_mutation_rate;
        auto ln_prob_misaligned = ln_probability_read_correctly_aligned(penalty, read, mu);
        return ln_prob_misaligned < options_.misalignment_parameters->min_ln_prob_correctly_aligned;
    } else {
        return false;
    }
}

CigarScanner::VariantObservation
CigarScanner::make_observation(const CandidateIterator first_match, const CandidateIterator last_match) const
{
    assert(first_match != last_match);
    const Candidate& candidate {*first_match};
    VariantObservation result {};
    result.variant = candidate.variant;
    result.total_depth = get_min_depth(candidate.variant, combined_read_coverage_tracker_);
    std::vector<Candidate> observations {first_match, last_match};
    std::sort(begin(observations), end(observations),
              [] (const Candidate& lhs, const Candidate& rhs) { return lhs.origin.get() < rhs.origin.get(); });
    for (auto observation_itr = begin(observations); observation_itr != end(observations);) {
        const auto& origin = observation_itr->origin;
        auto next_itr = std::find_if_not(next(observation_itr), end(observations),
                                         [&] (const Candidate& c) { return c.origin.get() == origin.get(); });
        const auto num_observations = static_cast<unsigned>(std::distance(observation_itr, next_itr));
        std::vector<unsigned> observed_base_qualities(num_observations);
        std::transform(observation_itr, next_itr, begin(observed_base_qualities),
                       [this] (const Candidate& c) noexcept { return sum_base_qualities(c); });
        std::vector<AlignedRead::MappingQuality> observed_mapping_qualities(num_observations);
        std::transform(observation_itr, next_itr, begin(observed_mapping_qualities),
                       [] (const Candidate& c) noexcept { return c.source.get().mapping_quality(); });
        const auto forward_strand_support = std::accumulate(observation_itr, next_itr, 0u,
                                                     [] (unsigned curr, const Candidate& c) noexcept {
                                                         if (is_forward_strand(c.source.get())) {
                                                             ++curr;
                                                         }
                                                         return curr;
                                                     });
        const auto edge_support = std::accumulate(observation_itr, next_itr, 0u,
                                                      [] (unsigned curr, const Candidate& c) noexcept {
                                                          if (begins_equal(c, c.source.get()) || ends_equal(c, c.source.get())) {
                                                              ++curr;
                                                          }
                                                          return curr;
                                                      });
        const auto depth = std::max(get_min_depth(candidate.variant, sample_read_coverage_tracker_.at(origin)), num_observations);
        const auto forward_depth = get_min_depth(candidate.variant, sample_forward_strand_coverage_tracker_.at(origin));
        result.sample_observations.push_back({origin, depth, forward_depth,
                                              std::move(observed_base_qualities),
                                              std::move(observed_mapping_qualities),
                                              forward_strand_support, edge_support});
        observation_itr = next_itr;
    }
    return result;
}

std::vector<Variant>
CigarScanner::get_novel_likely_misaligned_candidates(const std::vector<Variant>& current_candidates) const
{
    std::is_sorted(std::cbegin(likely_misaligned_candidates_), std::cend(likely_misaligned_candidates_));
    std::vector<Candidate> unique_misaligned_candidates {};
    unique_misaligned_candidates.reserve(likely_misaligned_candidates_.size());
    std::unique_copy(std::cbegin(likely_misaligned_candidates_), std::cend(likely_misaligned_candidates_),
                     std::back_inserter(unique_misaligned_candidates));
    std::vector<Variant> unique_misaligned_variants {};
    unique_misaligned_variants.reserve(unique_misaligned_candidates.size());
    std::transform(std::cbegin(unique_misaligned_candidates), std::cend(unique_misaligned_candidates),
                   std::back_inserter(unique_misaligned_variants),
                   [] (const Candidate& candidate) { return candidate.variant; });
    std::vector<Variant> result {};
    result.reserve(unique_misaligned_variants.size());
    assert(std::is_sorted(std::cbegin(current_candidates), std::cend(current_candidates)));
    std::set_difference(std::cbegin(unique_misaligned_variants), std::cend(unique_misaligned_variants),
                        std::cbegin(current_candidates), std::cend(current_candidates),
                        std::back_inserter(result));
    return result;
}

// non-member methods

namespace {

struct StrandSupportStats
{
    unsigned forward_support, forward_depth, reverse_support, reverse_depth;
};

auto sum(const std::vector<unsigned>& observed_qualities) noexcept
{
    return std::accumulate(std::cbegin(observed_qualities), std::cend(observed_qualities), 0);
}

void erase_below(std::vector<unsigned>& observed_qualities, const unsigned min)
{
    observed_qualities.erase(std::remove_if(std::begin(observed_qualities), std::end(observed_qualities),
                                            [=] (const auto q) { return q < min; }),
                             std::end(observed_qualities));
}

void partial_sort(std::vector<unsigned>& observed_qualities, const unsigned n)
{
    std::partial_sort(std::begin(observed_qualities), std::next(std::begin(observed_qualities), n),
                      std::end(observed_qualities), std::greater<> {});
}

double compute_strand_bias(const StrandSupportStats& strand_depths)
{
    return 1 - maths::fisher_exact_test(strand_depths.forward_support, strand_depths.forward_depth - strand_depths.forward_support,
                                        strand_depths.reverse_support, strand_depths.reverse_depth - strand_depths.reverse_support);
}

bool only_observed_on_one_strand(const StrandSupportStats& strand_depths) noexcept
{
    const auto support = strand_depths.forward_support + strand_depths.reverse_support;
    return support > 0 && (strand_depths.forward_support == 0 || strand_depths.reverse_support == 0);
}

bool is_likely_runthrough_artifact(const StrandSupportStats& strand_depths, std::vector<unsigned>& observed_qualities)
{
    const auto num_observations = strand_depths.forward_support + strand_depths.reverse_support;
    if (num_observations < 10 || !only_observed_on_one_strand(strand_depths)) return false;
    assert(!observed_qualities.empty());
    const auto median_bq = maths::median(observed_qualities);
    return median_bq < 15;
}

bool is_tandem_repeat(const Allele& allele, const unsigned max_period = 4)
{
    for (unsigned period {0}; period <= max_period; ++period) {
        if (utils::is_tandem_repeat(allele.sequence(), period)) return true;
    }
    return false;
}

bool is_good_germline(const Variant& variant, const unsigned depth, const unsigned forward_strand_depth,
                      const unsigned forward_strand_support, std::vector<unsigned> observed_qualities,
                      const unsigned copy_number = 2)
{
    const auto support = observed_qualities.size();
    if (depth < 4) {
        return support > 1 || sum(observed_qualities) >= 30 || is_deletion(variant);
    }
    const StrandSupportStats strand_depths {forward_strand_support,
                                            forward_strand_depth,
                                            static_cast<unsigned>(support) - forward_strand_support,
                                            depth - forward_strand_depth};
    const auto strand_bias = compute_strand_bias(strand_depths);
    if (support > 20 && strand_bias > 0.99 && only_observed_on_one_strand(strand_depths)) return false;
    if (is_snv(variant)) {
        if (is_likely_runthrough_artifact(strand_depths, observed_qualities)) return false;
        erase_below(observed_qualities, 20);
        if (depth <= 10) return observed_qualities.size() > 1;
        return observed_qualities.size() > 2 && static_cast<double>(observed_qualities.size()) / depth > (1. / (5 * copy_number));
    } else if (is_insertion(variant)) {
        if (support == 1 && alt_sequence_size(variant) > 10) return false;
        if (depth < 10) {
            return support > 1 || (alt_sequence_size(variant) > 3 && is_tandem_repeat(variant.alt_allele()));
        } else if (depth <= 30) {
            return support > 1;
        } else if (depth <= 60) {
            if (support == 1) return false;
            if (static_cast<double>(support) / depth > 0.3) return true;
            erase_below(observed_qualities, 25);
            if (observed_qualities.size() <= 1) return false;
            if (observed_qualities.size() > 2) return true;
            partial_sort(observed_qualities, 2);
            return static_cast<double>(observed_qualities[0]) / alt_sequence_size(variant) > 20;
        } else {
            if (support == 1) return false;
            if (static_cast<double>(support) / depth > 0.35) return true;
            erase_below(observed_qualities, 20);
            if (observed_qualities.size() <= 1) return false;
            if (observed_qualities.size() > 3) return true;
            return static_cast<double>(observed_qualities[0]) / alt_sequence_size(variant) > 20;
        }
    } else {
        // deletion or mnv
        if (region_size(variant) < 10) {
            return support > 1 && static_cast<double>(support) / depth > (1. / (10 * copy_number));
        } else {
            return static_cast<double>(support) / (depth - std::sqrt(depth)) > (1. / (5 * copy_number));
        }
    }
}

struct UnknownExpectedVAFStats
{
    double min_vaf, min_probability;
    AlignedRead::BaseQuality min_bq = 15;
};

auto beta_sf(unsigned a, unsigned b, double x)
{
    // Haldane's prior but make sure is proper
    return maths::beta_sf(static_cast<double>(std::max(a, 1u)), static_cast<double>(std::max(b, 1u)), x);
}

bool is_good_somatic(const Variant& variant, const unsigned depth, const unsigned forward_strand_depth,
                     const unsigned forward_strand_support, const unsigned num_edge_observations,
                     std::vector<unsigned> observed_qualities, const UnknownExpectedVAFStats vaf_def)
{
    assert(depth > 0);
    const auto support = observed_qualities.size();
    const StrandSupportStats strand_depths {forward_strand_support,
                                            forward_strand_depth,
                                            static_cast<unsigned>(support) - forward_strand_support,
                                            depth - forward_strand_depth};
    const auto strand_bias = compute_strand_bias(strand_depths);
    const auto raw_vaf = static_cast<double>(support) / depth;
    if (support > 10 && strand_bias > 0.99) {
        if (only_observed_on_one_strand(strand_depths)) return false;
        if (strand_bias > 0.99999999 && raw_vaf < 0.5) return false;
    }
    if (is_snv(variant)) {
        if (is_likely_runthrough_artifact(strand_depths, observed_qualities)) return false;
        erase_below(observed_qualities, vaf_def.min_bq);
        if (observed_qualities.size() <= num_edge_observations) return false;
        const auto good_support = observed_qualities.size() - num_edge_observations;
        const auto probability_vaf_greater_than_min_vaf = beta_sf(good_support, depth - good_support, vaf_def.min_vaf);
        return good_support > 1
            && probability_vaf_greater_than_min_vaf >= vaf_def.min_probability
            && num_edge_observations < support;
    } else if (is_insertion(variant)) {
        if (support == 1 && alt_sequence_size(variant) > 8) return false;
        erase_below(observed_qualities, vaf_def.min_bq);
        const auto good_support = observed_qualities.size();
        if (good_support > 1 && alt_sequence_size(variant) > 10) return true;
        const auto probability_vaf_greater_than_min_vaf = beta_sf(good_support, depth - good_support, vaf_def.min_vaf);
        return good_support > 1 && probability_vaf_greater_than_min_vaf >= vaf_def.min_probability;
    } else {
        // deletion or mnv
        const auto probability_vaf_greater_than_min_vaf = beta_sf(support, depth - support, vaf_def.min_vaf);
        return support > 1 && probability_vaf_greater_than_min_vaf >= vaf_def.min_probability;
    }
}

bool is_good_germline(const Variant& v, const CigarScanner::VariantObservation::SampleObservationStats& observation,
                      unsigned ploidy = 2)
{
    return is_good_germline(v, observation.depth, observation.forward_strand_depth,
                            observation.forward_strand_support, observation.observed_base_qualities,
                            ploidy);
}

bool any_good_germline_samples(const CigarScanner::VariantObservation& candidate, unsigned ploidy = 2)
{
    return std::any_of(std::cbegin(candidate.sample_observations), std::cend(candidate.sample_observations),
                       [&] (const auto& observation) { return is_good_germline(candidate.variant, observation, ploidy); });
}

auto count_forward_strand_depth(const CigarScanner::VariantObservation& candidate)
{
    return std::accumulate(std::cbegin(candidate.sample_observations), std::cend(candidate.sample_observations), 0u,
                           [&] (auto curr, const auto& observation) { return curr + observation.forward_strand_depth; });
}

auto count_forward_strand_support(const CigarScanner::VariantObservation& candidate)
{
    return std::accumulate(std::cbegin(candidate.sample_observations), std::cend(candidate.sample_observations), 0u,
                           [&] (auto curr, const auto& observation) { return curr + observation.forward_strand_support; });
}

auto count_edge_support(const CigarScanner::VariantObservation& candidate)
{
    return std::accumulate(std::cbegin(candidate.sample_observations), std::cend(candidate.sample_observations), 0u,
                           [&] (auto curr, const auto& observation) { return curr + observation.edge_support; });
}

auto concat_observed_base_qualities(const CigarScanner::VariantObservation& candidate)
{
    std::size_t num_base_qualities {0};
    for (const auto& observation : candidate.sample_observations) {
        num_base_qualities += observation.observed_base_qualities.size();
    }
    std::vector<unsigned> result {};
    result.reserve(num_base_qualities);
    for (const auto& observation : candidate.sample_observations) {
        utils::append(observation.observed_base_qualities, result);
    }
    return result;
}

bool is_good_germline_pooled(const CigarScanner::VariantObservation& candidate)
{
    return is_good_germline(candidate.variant, candidate.total_depth, count_forward_strand_depth(candidate),
                            count_forward_strand_support(candidate), concat_observed_base_qualities(candidate));
}

bool is_good_somatic(const Variant& v, const CigarScanner::VariantObservation::SampleObservationStats& observation,
                     UnknownExpectedVAFStats vaf_def)
{
    return is_good_somatic(v, observation.depth, observation.forward_strand_depth, observation.forward_strand_support,
                           observation.edge_support, observation.observed_base_qualities, vaf_def);
}

bool is_good_somatic_pooled(const CigarScanner::VariantObservation& candidate, const UnknownExpectedVAFStats& vaf_def)
{
    return is_good_somatic(candidate.variant, candidate.total_depth, count_forward_strand_depth(candidate),
                           count_forward_strand_support(candidate), count_edge_support(candidate),
                           concat_observed_base_qualities(candidate), vaf_def);
}

bool is_good_pacbio(const Variant& variant, const unsigned depth, const unsigned forward_strand_depth,
                    const unsigned forward_strand_support, std::vector<unsigned> observed_qualities)
{
    const auto support = observed_qualities.size();
    const auto vaf = static_cast<double>(support) / depth;
    if (is_snv(variant)) {
        return support > 1 && vaf > 0.1;
    } else if (is_insertion(variant)) {
        if (alt_sequence_size(variant) > 500) {
            return true;
        } else if (alt_sequence_size(variant) > 200) {
            return vaf > 0.02;
        } else if (alt_sequence_size(variant) > 20) {
            return vaf > 0.05;
        }
        if (support < 2) return false;
        if (alt_sequence_size(variant) <= 2) {
            return vaf > 0.2;
        } else if (alt_sequence_size(variant) < 4) {
            return vaf > 0.1;
        } else {
            return vaf > 0.05;
        }
    } else { // deletion or mnv
        if (region_size(variant) > 50) {
            return vaf > 0.1;
        }
        if (support < 2) return false;
        if (region_size(variant) <= 2) {
            return vaf > 0.2;
        } else if (region_size(variant) < 4) {
            return vaf > 0.1;
        } else {
            return vaf > 0.05;
        }
    }
}

bool is_good_pacbio(const Variant& v, const CigarScanner::VariantObservation::SampleObservationStats& observation)
{
    return is_good_pacbio(v, observation.depth, observation.forward_strand_depth,
                          observation.forward_strand_support, observation.observed_base_qualities);
}

bool any_good_pacbio_samples(const CigarScanner::VariantObservation& candidate)
{
    return std::any_of(std::cbegin(candidate.sample_observations), std::cend(candidate.sample_observations),
                       [&] (const auto& observation) { return is_good_pacbio(candidate.variant, observation); });
}

bool is_good_pacbio_pooled(const CigarScanner::VariantObservation& candidate)
{
    return is_good_pacbio(candidate.variant, candidate.total_depth, count_forward_strand_depth(candidate),
                          count_forward_strand_support(candidate), concat_observed_base_qualities(candidate));
}

} // namespace

bool KnownCopyNumberInclusionPredicate::operator()(const CigarScanner::VariantObservation& candidate)
{
    return any_good_germline_samples(candidate, copy_number_) || (candidate.sample_observations.size() > 1 && is_good_germline_pooled(candidate));
}

bool PacBioInclusionPredicate::operator()(const CigarScanner::VariantObservation& candidate)
{
    return any_good_pacbio_samples(candidate) || (candidate.sample_observations.size() > 1 && is_good_pacbio_pooled(candidate));
}

UnknownCopyNumberInclusionPredicate::UnknownCopyNumberInclusionPredicate(double min_vaf, double min_probability, AlignedRead::BaseQuality min_bq)
: normal_ {}
, min_vaf_ {min_vaf}
, min_probability_ {min_probability}
, min_bq_ {min_bq}
{}

UnknownCopyNumberInclusionPredicate::UnknownCopyNumberInclusionPredicate(SampleName normal, double min_vaf, double min_probability, AlignedRead::BaseQuality min_bq)
: normal_ {std::move(normal)}
, min_vaf_ {min_vaf}
, min_probability_ {min_probability}
, min_bq_ {min_bq}
{}

bool UnknownCopyNumberInclusionPredicate::operator()(const CigarScanner::VariantObservation& candidate)
{
    const UnknownExpectedVAFStats vaf_def {min_vaf_, min_probability_, min_bq_};
    return std::any_of(std::cbegin(candidate.sample_observations), std::cend(candidate.sample_observations),
                       [&] (const auto& observation) {
                           if (normal_ && observation.sample.get() == *normal_) {
                               return is_good_germline(candidate.variant, observation);
                           } else {
                               return is_good_somatic(candidate.variant, observation, vaf_def);
                           }
                       });
}

bool is_good_cell(const Variant& v, const CigarScanner::VariantObservation::SampleObservationStats& observation)
{
    const UnknownExpectedVAFStats vaf_def {0.2, 0.5};
    return is_good_somatic(v, observation, vaf_def);
}

bool any_good_cell_samples(const CigarScanner::VariantObservation& candidate)
{
    return std::any_of(std::cbegin(candidate.sample_observations), std::cend(candidate.sample_observations),
                       [&] (const auto& observation) { return is_good_cell(candidate.variant, observation); });
}

bool is_good_cell_pooled(const CigarScanner::VariantObservation& candidate)
{
    const UnknownExpectedVAFStats vaf_def {0.2, 0.8};
    return is_good_somatic_pooled(candidate, vaf_def);
}

bool CellInclusionPredicate::operator()(const CigarScanner::VariantObservation& candidate)
{
    return any_good_cell_samples(candidate) || (candidate.sample_observations.size() > 1 && is_good_cell_pooled(candidate));
}

namespace {

auto count_observations(const CigarScanner::VariantObservation& candidate)
{
    return std::accumulate(std::cbegin(candidate.sample_observations), std::cend(candidate.sample_observations), std::size_t {0},
                           [] (auto curr, const auto& sample) { return curr + sample.observed_base_qualities.size(); });
}

} // namespace

bool SimpleThresholdInclusionPredicate::operator()(const CigarScanner::VariantObservation& candidate) noexcept
{
    return count_observations(candidate) >= min_observations_;
}

bool TolerantMatchPredicate::operator()(const Variant& lhs, const Variant& rhs) noexcept
{
    if (!are_same_type(lhs, rhs) || is_snv(lhs) || is_mnv(lhs)) {
        return lhs == rhs;
    }
    if (is_insertion(lhs)) {
        if (!is_same_region(lhs, rhs)) return false;
        if (alt_sequence_size(lhs) == alt_sequence_size(rhs)) {
            const auto& lhs_alt = alt_sequence(lhs);
            const auto& rhs_alt = alt_sequence(rhs);
            return std::count(std::cbegin(lhs_alt), std::cend(lhs_alt), 'N')
                   == std::count(std::cbegin(rhs_alt), std::cend(rhs_alt), 'N');
        } else {
            if (utils::is_homopolymer(alt_sequence(lhs)) && utils::is_homopolymer(alt_sequence(rhs))) {
                return alt_sequence(lhs).front() == alt_sequence(rhs).front();
            } else {
                return false;
            }
        }
    }
    return overlaps(lhs, rhs);
}

} // coretools
} // namespace octopus
