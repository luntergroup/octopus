// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "dynamic_cigar_scanner.hpp"

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
#include "logging/logging.hpp"

#include "utils/maths.hpp"

namespace octopus { namespace coretools {

std::unique_ptr<VariantGenerator> DynamicCigarScanner::do_clone() const
{
    return std::make_unique<DynamicCigarScanner>(*this);
}

DynamicCigarScanner::DynamicCigarScanner(const ReferenceGenome& reference, Options options)
: reference_ {reference}
, options_ {options}
, buffer_ {}
, candidates_ {}
, likely_misaligned_candidates_ {}
, max_seen_candidate_size_ {}
, read_coverage_tracker_ {}
, misaligned_tracker_ {}
{
    buffer_.reserve(100);
}

bool DynamicCigarScanner::do_requires_reads() const noexcept
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

void DynamicCigarScanner::add_read(const AlignedRead& read)
{
    using std::cbegin; using std::next; using std::move;
    using Flag = CigarOperation::Flag;
    
    const auto& read_contig   = contig_name(read);
    const auto& read_sequence = read.sequence();
    auto sequence_iter     = cbegin(read_sequence);
    auto base_quality_iter = cbegin(read.base_qualities());
    auto ref_index = mapped_begin(read);
    std::size_t read_index {0};
    GenomicRegion region;
    double misalignment_penalty {0};
    for (const auto& cigar_operation : read.cigar()) {
        const auto op_size = cigar_operation.size();
        switch (cigar_operation.flag()) {
        case Flag::alignmentMatch:
            misalignment_penalty += add_snvs_in_match_range(GenomicRegion {read_contig, ref_index, ref_index + op_size},
                                                            next(sequence_iter, read_index),
                                                            next(sequence_iter, read_index + op_size),
                                                            next(base_quality_iter, read_index));
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
            add_candidate(region,
                          reference_.get().fetch_sequence(region),
                          copy(read_sequence, read_index, op_size),
                          next(base_quality_iter, read_index));
            read_index += op_size;
            ref_index  += op_size;
            misalignment_penalty += op_size * options_.misalignment_parameters.snv_penalty;
            break;
        }
        case Flag::insertion:
        {
            add_candidate(GenomicRegion {read_contig, ref_index, ref_index},
                          "",
                          copy(read_sequence, read_index, op_size),
                          next(base_quality_iter, read_index));
            read_index += op_size;
            misalignment_penalty += options_.misalignment_parameters.indel_penalty;
            break;
        }
        case Flag::deletion:
        {
            region = GenomicRegion {read_contig, ref_index, ref_index + op_size};
            add_candidate(move(region),
                          reference_.get().fetch_sequence(region),
                          "",
                          next(base_quality_iter, read_index));
            ref_index += op_size;
            misalignment_penalty += options_.misalignment_parameters.indel_penalty;
            break;
        }
        case Flag::softClipped:
        {
            read_index += op_size;
            ref_index  += op_size;
            if (op_size > options_.misalignment_parameters.max_unpenalised_clip_size) {
                misalignment_penalty += options_.misalignment_parameters.clip_penalty;
            }
            break;
        }
        case Flag::hardClipped:
        {
            if (op_size > options_.misalignment_parameters.max_unpenalised_clip_size) {
                misalignment_penalty += options_.misalignment_parameters.clip_penalty;
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
        read_coverage_tracker_.add(clipped_mapped_region(read));
    } else {
        read_coverage_tracker_.add(read);
    }
    const auto mu = options_.misalignment_parameters.max_expected_mutation_rate;
    const auto ln_prob_misaligned = ln_probability_read_correctly_aligned(misalignment_penalty, read, mu);
    const auto min_ln_prob_misaligned = options_.misalignment_parameters.min_ln_prob_correctly_aligned;
    if (ln_prob_misaligned >= min_ln_prob_misaligned) {
        utils::append(std::move(buffer_), candidates_);
    } else {
        utils::append(std::move(buffer_), likely_misaligned_candidates_);
        misaligned_tracker_.add(clipped_mapped_region(read));
    }
}

void DynamicCigarScanner::do_add_read(const SampleName& sample, const AlignedRead& read)
{
    add_read(read);
}

void DynamicCigarScanner::do_add_reads(const SampleName& sample, VectorIterator first, VectorIterator last)
{
    std::for_each(first, last, [this] (const AlignedRead& read) { add_read(read); });
}

void DynamicCigarScanner::do_add_reads(const SampleName& sample, FlatSetIterator first, FlatSetIterator last)
{
    std::for_each(first, last, [this] (const AlignedRead& read) { add_read(read); });
}

unsigned get_min_depth(const Variant& v, const CoverageTracker<GenomicRegion>& tracker)
{
    if (is_insertion(v)) {
        const auto& region = mapped_region(v);
        if (region.begin() > 0) {
            return tracker.min_coverage(expand(region, 1, 1));
        } else {
            return tracker.min_coverage(expand_rhs(region, 1));
        }
    } else {
        return tracker.min_coverage(mapped_region(v));
    }
}

struct RegionVariants : public Mappable<RegionVariants>
{
    RegionVariants(GenomicRegion region) : region {std::move(region)} {}
    GenomicRegion region;
    std::deque<Variant> variants;
    const GenomicRegion& mapped_region() const noexcept { return region; }
};

auto init_region_variants(const std::vector<GenomicRegion>& regions)
{
    std::vector<RegionVariants> result {};
    result.reserve(regions.size());
    std::transform(std::cbegin(regions), std::cend(regions), std::back_inserter(result),
                   [] (const auto& region) { return RegionVariants {region}; });
    return result;
}

boost::optional<RegionVariants&> find_contained(std::vector<RegionVariants>& regions, const Variant& variant)
{
    if (regions.empty()) return boost::none;
    const auto itr = std::find_if(std::begin(regions), std::end(regions),
                                  [&] (const auto& region) { return contains(region, variant); });
    if (itr != std::end(regions)) {
        assert(contains(*itr, variant));
        return *itr;
    } else {
        return boost::none;
    }
}

std::vector<Variant> DynamicCigarScanner::do_generate_variants(const GenomicRegion& region)
{
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::next; using std::distance;
    
    std::sort(begin(candidates_), end(candidates_));
    auto overlapped = overlap_range(candidates_, region, max_seen_candidate_size_);
    
    std::vector<Variant> result {};
    
    if (empty(overlapped)) return result;
    
    result.reserve(size(overlapped, BidirectionallySortedTag {})); // maximum possible
    
    const auto repeat_regions = get_repeat_regions(region);
    auto repeat_region_variants = init_region_variants(repeat_regions);
    
    while (!overlapped.empty()) {
        const Candidate& candidate {overlapped.front()};
        const auto next_candidate = std::find_if_not(next(cbegin(overlapped)), cend(overlapped),
                                                     [this, &candidate] (const Candidate& c) {
                                                         return options_.match(c.variant, candidate.variant);
                                                     });
        const auto num_observations = static_cast<unsigned>(distance(cbegin(overlapped), next_candidate));
        const auto min_depth = get_min_depth(candidate.variant, read_coverage_tracker_);
        thread_local std::vector<unsigned> observed_qualities {};
        observed_qualities.resize(num_observations);
        std::transform(cbegin(overlapped), next_candidate, begin(observed_qualities),
                       [this] (const Candidate& c) noexcept { return sum_base_qualities(c); });
        if (options_.include(candidate.variant, min_depth, observed_qualities)) {
            if (num_observations > 1) {
                auto unique_iter = cbegin(overlapped);
                while (unique_iter != next_candidate) {
                    auto repeat_region = find_contained(repeat_region_variants, unique_iter->variant);
                    if (repeat_region) {
                        repeat_region->variants.push_back(unique_iter->variant);
                    } else {
                        result.push_back(unique_iter->variant);
                    }
                    unique_iter = std::find_if_not(next(unique_iter), next_candidate,
                                                   [unique_iter] (const Candidate& c) {
                                                       return c.variant == unique_iter->variant;
                                                   });
                }
            } else {
                auto repeat_region = find_contained(repeat_region_variants, candidate.variant);
                if (repeat_region) {
                    repeat_region->variants.push_back(candidate.variant);
                } else {
                    result.push_back(candidate.variant);
                }
            }
        }
        overlapped.advance_begin(num_observations);
    }
    
    std::sort(std::begin(likely_misaligned_candidates_), std::end(likely_misaligned_candidates_));
    std::vector<Candidate> unique_misaligned_candidates {};
    unique_misaligned_candidates.reserve(likely_misaligned_candidates_.size());
    std::unique_copy(std::cbegin(likely_misaligned_candidates_), std::cend(likely_misaligned_candidates_),
                     std::back_inserter(unique_misaligned_candidates));
    std::vector<Variant> unique_misaligned_variants {};
    unique_misaligned_variants.reserve(unique_misaligned_candidates.size());
    std::transform(std::cbegin(unique_misaligned_candidates), std::cend(unique_misaligned_candidates),
                   std::back_inserter(unique_misaligned_variants),
                   [] (const Candidate& candidate) { return candidate.variant; });
    std::vector<Variant> novel_unique_misaligned_variants {};
    novel_unique_misaligned_variants.reserve(unique_misaligned_variants.size());
    std::set_difference(std::cbegin(unique_misaligned_variants), std::cend(unique_misaligned_variants),
                        std::cbegin(result), std::cend(result), std::back_inserter(novel_unique_misaligned_variants));
    if (!novel_unique_misaligned_variants.empty() && debug_log_) {
        stream(*debug_log_) << "DynamicCigarScanner: ignoring " << count_overlapped(novel_unique_misaligned_variants, region)
                            << " unique candidates in " << region;
    }
    for (const auto& candidate : novel_unique_misaligned_variants) {
        auto repeat_region = find_contained(repeat_region_variants, candidate);
        if (repeat_region) {
            repeat_region->variants.push_back(candidate);
        }
    }
    
    std::vector<Variant> good_repeat_region_variants {};
    for (auto& region : repeat_region_variants) {
        const auto density = region.variants.size() / size(region.region);
        if (density < options_.max_repeat_region_density) {
            utils::append(std::move(region.variants), good_repeat_region_variants);
        } else {
            if (debug_log_) {
                stream(*debug_log_) << "DynamicCigarScanner: masking " << region.variants.size()
                                    << " candidates in repetitive region " << region.region;
            }
        }
    }
    auto itr = utils::append(std::move(good_repeat_region_variants), result);
    std::inplace_merge(std::begin(result), itr, std::end(result));
    
    return result;
}

void DynamicCigarScanner::do_clear() noexcept
{
    buffer_.clear();
    buffer_.shrink_to_fit();
    candidates_.clear();
    candidates_.shrink_to_fit();
    likely_misaligned_candidates_.clear();
    likely_misaligned_candidates_.shrink_to_fit();
    read_coverage_tracker_.clear();
    misaligned_tracker_.clear();
    max_seen_candidate_size_ = 0;
}

std::string DynamicCigarScanner::name() const
{
    return "DynamicCigarScanner";
}

// private methods

double DynamicCigarScanner::add_snvs_in_match_range(const GenomicRegion& region,
                                                    const SequenceIterator first_base, const SequenceIterator last_base,
                                                    AlignedRead::BaseQualityVector::const_iterator first_base_quality)
{
    using boost::make_zip_iterator; using std::for_each; using std::cbegin; using std::cend;
    using Tuple = boost::tuple<char, char>;
    const NucleotideSequence ref_segment {reference_.get().fetch_sequence(region)};
    const auto& contig = region.contig_name();
    auto ref_index = mapped_begin(region);
    double misalignment_penalty {0};
    for_each(make_zip_iterator(boost::make_tuple(cbegin(ref_segment), first_base)),
             make_zip_iterator(boost::make_tuple(cend(ref_segment), last_base)),
             [this, &contig, &ref_index, &first_base_quality, &misalignment_penalty] (const Tuple& t) {
                 const char ref_base  {t.get<0>()}, read_base {t.get<1>()};
                 if (ref_base != read_base && ref_base != 'N' && read_base != 'N') {
                     add_candidate(GenomicRegion {contig, ref_index, ref_index + 1},
                                   ref_base, read_base, first_base_quality);
                     if (*first_base_quality >= options_.misalignment_parameters.snv_threshold) {
                         misalignment_penalty += options_.misalignment_parameters.snv_penalty;
                     }
                 }
                 ++ref_index;
                 ++first_base_quality;
             });
    return misalignment_penalty;
}

unsigned DynamicCigarScanner::sum_base_qualities(const Candidate& candidate) const noexcept
{
    return std::accumulate(candidate.first_base_quality_iter,
                           std::next(candidate.first_base_quality_iter, alt_sequence_size(candidate.variant)),
                           0u);
}

std::vector<GenomicRegion> DynamicCigarScanner::get_repeat_regions(const GenomicRegion& region) const
{
    if (options_.repeat_region_generator) {
        return (*options_.repeat_region_generator)(reference_, region);
    } else {
        return {};
    }
}

// non-member methods

namespace {

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
    std::partial_sort(std::begin(observed_qualities), std::next(std::begin(observed_qualities), 2),
                      std::end(observed_qualities), std::greater<> {});
}

} // namespace

bool DefaultInclusionPredicate::operator()(const Variant& v, const unsigned depth, std::vector<unsigned>& observed_qualities)
{
    const auto num_observations = observed_qualities.size();
    if (depth < 4) {
        return num_observations > 1 || sum(observed_qualities) >= 20 || is_deletion(v);
    }
    if (is_snp(v)) {
        const auto base_quality_sum = sum(observed_qualities);
        if (depth <= 60) {
            if (num_observations < 2) return false;
            if (base_quality_sum > 100) return true;
            erase_below(observed_qualities, 5);
            if (observed_qualities.size() < 2) return false;
            if (static_cast<double>(observed_qualities.size()) / depth > 0.2) return true;
            partial_sort(observed_qualities, 2);
            return observed_qualities[0] >= 20 && observed_qualities[1] >= 20;
        } else {
            if (num_observations < 3) return false;
            if (base_quality_sum > 150) return true;
            erase_below(observed_qualities, 10);
            if (observed_qualities.size() < 3) return false;
            if (static_cast<double>(observed_qualities.size()) / depth > 0.2) return true;
            partial_sort(observed_qualities, 3);
            return observed_qualities[0] >= 30 && observed_qualities[1] >= 25 && observed_qualities[2] >= 20;
        }
    } else if (is_insertion(v)) {
        if (num_observations == 1 && alt_sequence_size(v) > 8) return false;
        if (depth <= 15) {
            return num_observations > 1;
        } else if (depth <= 30) {
            if (static_cast<double>(num_observations) / depth > 0.45) return true;
            erase_below(observed_qualities, 20);
            return num_observations > 1;
        } else if (depth <= 60) {
            if (num_observations == 1) return false;
            if (static_cast<double>(num_observations) / depth > 0.4) return true;
            erase_below(observed_qualities, 25);
            if (observed_qualities.size() <= 1) return false;
            if (observed_qualities.size() > 2) return true;
            partial_sort(observed_qualities, 2);
            return static_cast<double>(observed_qualities[0]) / alt_sequence_size(v) > 20;
        } else {
            if (num_observations == 1) return false;
            if (static_cast<double>(num_observations) / depth > 0.35) return true;
            erase_below(observed_qualities, 20);
            if (observed_qualities.size() <= 1) return false;
            if (observed_qualities.size() > 3) return true;
            return static_cast<double>(observed_qualities[0]) / alt_sequence_size(v) > 20;
        }
    } else {
        // deletion or mnv
        if (region_size(v) < 10) {
            return num_observations > 1 && static_cast<double>(num_observations) / depth > 0.05;
        } else {
            return num_observations > 1 && num_observations >= 5;
        }
    }
}

bool DefaultMatchPredicate::operator()(const Variant& lhs, const Variant& rhs) noexcept
{
    if (!are_same_type(lhs, rhs) || is_snp(lhs) || is_mnv(lhs)) {
        return lhs == rhs;
    }
    if (is_insertion(lhs) && alt_sequence_size(lhs) == alt_sequence_size(rhs)) {
        const auto& lhs_alt = alt_sequence(lhs);
        const auto& rhs_alt = alt_sequence(rhs);
        return std::count(std::cbegin(lhs_alt), std::cend(lhs_alt), 'N')
               == std::count(std::cbegin(rhs_alt), std::cend(rhs_alt), 'N');
    }
    return overlaps(lhs, rhs);
}

} // coretools
} // namespace octopus
