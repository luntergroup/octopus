// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "assembler_active_region_generator.hpp"

#include <iterator>
#include <algorithm>
#include <cmath>
#include <cassert>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>

#include "basics/genomic_region.hpp"
#include "basics/cigar_string.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/maths.hpp"
#include "utils/append.hpp"
#include "utils/free_memory.hpp"

#include <iostream>

namespace octopus { namespace coretools {

AssemblerActiveRegionGenerator::AssemblerActiveRegionGenerator(const ReferenceGenome& reference)
: reference_ {reference}
, coverage_tracker_ {}
, interesting_read_coverages_ {}
, clipped_coverage_tracker_ {}
{}

namespace {

template <typename Container>
bool contains(const Container& values, const typename Container::value_type& value)
{
    return std::find(std::cbegin(values), std::cend(values), value) != std::cend(values);
}

} // namespace

AssemblerActiveRegionGenerator::AssemblerActiveRegionGenerator(const ReferenceGenome& reference, Options options)
: reference_ {reference}
, coverage_tracker_ {}
, interesting_read_coverages_ {}
, clipped_coverage_tracker_ {}
{
    snvs_interesting_ = contains(options.trigger_types, Options::TriggerType::snv);
    indels_interesting_ = contains(options.trigger_types, Options::TriggerType::indel);
    structual_interesting_ = contains(options.trigger_types, Options::TriggerType::structual);
	clustered_interesting_ = contains(options.trigger_types, Options::TriggerType::clustered);
    trigger_quality_ = options.trigger_quality;
    trigger_clip_size_ = options.trigger_clip_size;
    min_expected_mutation_frequency_ = options.min_expected_mutation_frequency;
    read_profile_ = options.read_profile;
}

void AssemblerActiveRegionGenerator::add(const SampleName& sample, const AlignedRead& read)
{
    coverage_tracker_[sample].add(read);
    if (is_interesting(read)) {
        interesting_read_coverages_[sample].add(read);
    }
    if (structual_interesting_) {
        if (is_soft_clipped(read)) {
            clipped_coverage_tracker_[sample].add(clipped_mapped_region(read));
        } else {
            clipped_coverage_tracker_[sample].add(read);
        }
    }
}

void AssemblerActiveRegionGenerator::add(const SampleName& sample, const AlignedTemplate& reads)
{
    for (const auto& read : reads) add(sample, read);
}

template <typename Container>
void log_each(Container& values)
{
    for (auto& value : values) value = std::log(value);
}

template <typename Container>
auto to_logs(const Container& values)
{
    auto result = values;
    log_each(result);
    return result;
}

struct HiddenStateParameters
{
    double mean, stdev;
};

template <typename RealType>
RealType log_normal_pdf(const RealType x, const RealType mean, const RealType stdev)
{
    constexpr RealType c {0.9189385332046727417803297364056176398613974736377834};
    return std::log(stdev) - c - std::pow(x - mean, 2) / (2 * std::pow(stdev, 2));
}

double emmision_probability(const unsigned x, const HiddenStateParameters params)
{
    return log_normal_pdf(static_cast<double>(x), params.mean, params.stdev);
}

auto compute_base_deletion_probabilities(const std::vector<unsigned>& coverages,
                                         const double mean_coverage,
                                         const double stdev_coverage)
{
    constexpr std::size_t num_states {3};
    constexpr std::array<double, num_states> begin_probabilities {1 - 2e-10, 1e-10, 1e-10};
    constexpr std::array<std::array<double, num_states + 1>, num_states> transition_probabilities {{
                                                                                                   {1.0 - (1e-2 + 1e-3 + 1e-20), 1e-2, 1e-3, 1e-20}, // normal
                                                                                                   {1.0 - (1e-7 + 1e-15 + 1e-20), 1e-7, 1e-15, 1e-20}, // insertion
                                                                                                   {1e-5, 1e-15, 1.0 - (1e-5 + 1e-15 + 1e-20), 1e-20}, // deletion
                                                                                                   }};
    constexpr std::array<double, num_states> end_probabilities {1 - 2e-10, 1e-10, 1e-10};
    std::array<HiddenStateParameters, num_states> state_params {{{mean_coverage, stdev_coverage},
                                                                {mean_coverage + 3 * stdev_coverage, stdev_coverage},
                                                                {std::max(mean_coverage - 2 * stdev_coverage, 0.0), stdev_coverage / 2}}};
    const auto num_observations = coverages.size();
    // convert to log space
    const auto a_0 = to_logs(begin_probabilities);
    const auto a_e = to_logs(end_probabilities);
    auto a = transition_probabilities;
    for (auto& c : a) log_each(c);
    using State = std::array<double, num_states>;
    using maths::log_sum_exp;
    // Forward part
    std::vector<State> fwd(num_observations);
    fwd[0] = a_0;
    for (std::size_t s {0}; s < num_states; ++s) {
        fwd[0][s] += emmision_probability(coverages[0], state_params[s]);
    }
    for (std::size_t i {1}; i < num_observations; ++i) {
        for (std::size_t s {0}; s < num_states; ++s) {
            State tmp {};
            for (std::size_t t {0}; t < num_states; ++t) {
                tmp[t] = fwd[i - 1][t] + a[t][s];
            }
            fwd[i][s] = log_sum_exp(tmp) + emmision_probability(coverages[i], state_params[s]);
        }
    }
    auto end = fwd.back();
    std::transform(std::cbegin(end), std::cend(end), std::cbegin(a_e), std::begin(end),
                   [] (auto p, auto l) { return p + l; });
    const auto p_fwd = log_sum_exp(end);
    // Backward part
    std::vector<State> bkw(num_observations);
    bkw.back() = a_e;
    for (int i = num_observations - 2; i >= 0; --i) {
        for (std::size_t s {0}; s < num_states; ++s) {
            State tmp {};
            for (std::size_t t {0}; t < num_states; ++t) {
                tmp[t] = bkw[i + 1][t] + emmision_probability(coverages[i + 1], state_params[t]) + a[s][t];
            }
            bkw[i][s] = log_sum_exp(tmp);
        }
    }
    // Merge
    std::vector<State> posteriors(num_observations);
    std::transform(std::cbegin(fwd), std::cend(fwd), std::cbegin(bkw), std::begin(posteriors),
                   [=] (const auto f, const auto b) {
                       State result {};
                       for (std::size_t s {0}; s < num_states; ++s) {
                           result[s] = std::exp(f[s] + b[s] - p_fwd);
                       }
                       return result;
                   });
    std::vector<double> result(num_observations);
    std::transform(std::cbegin(posteriors), std::cend(posteriors), std::begin(result),
                   [] (const auto& s) { return s[2]; });
    return result;
}

template <typename Container>
auto expand_each(const Container& regions, const GenomicRegion::Distance n)
{
    std::vector<GenomicRegion> result {};
    result.reserve(regions.size());
    std::transform(std::cbegin(regions), std::cend(regions), std::back_inserter(result),
                   [n] (const auto& region) { return expand(region, n); });
    return result;
}

auto get_deletion_hotspots(const GenomicRegion& region, const CoverageTracker<GenomicRegion>& tracker,
                           const boost::optional<const ReadSetProfile&>& read_profile = boost::none)
{
    const auto coverages = tracker.get(region);
    const auto mean_coverage = read_profile ? read_profile->depth_stats.combined.genome.all.mean : tracker.mean(region);
    const auto stdev_coverage = read_profile ? read_profile->depth_stats.combined.genome.all.stdev : tracker.stdev(region);
    const auto deletion_base_probs = compute_base_deletion_probabilities(coverages, mean_coverage, stdev_coverage);
    std::vector<bool> deletion_bases(deletion_base_probs.size());
    std::transform(std::cbegin(deletion_base_probs), std::cend(deletion_base_probs), std::begin(deletion_bases),
                   [] (const auto p) {  return p > 0.5; });
    return extract_covered_regions(expand_each(select_regions(region, deletion_bases), 50));
}

auto get_interesting_hotspots(const GenomicRegion& region,
                              const std::vector<unsigned>& interesting_coverages,
                              const std::vector<unsigned>& coverages,
                              const double min_expected_mutation_frequency)
{
    std::vector<bool> interesting_bases(interesting_coverages.size());
    std::transform(std::cbegin(coverages), std::cend(coverages), std::cbegin(interesting_coverages),
                   std::begin(interesting_bases),
                   [=] (const auto coverage, const auto interesting_coverage) {
                       if (coverage < 10) {
                           return interesting_coverage > 0;
                       } else if (coverage <= 30) {
                           return interesting_coverage > 1;
                       } else {
                           return static_cast<double>(interesting_coverage) / coverage >= min_expected_mutation_frequency;
                       }
                   });
    return select_regions(region, interesting_bases);
}

auto get_interesting_hotspots(const GenomicRegion& region,
                              const CoverageTracker<GenomicRegion>& interesting_read_tracker,
                              const CoverageTracker<GenomicRegion>& tracker,
                              const double min_expected_mutation_frequency)
{
    const auto interesting_coverages = interesting_read_tracker.get(region);
    const auto coverages = tracker.get(region);
    return get_interesting_hotspots(region, interesting_coverages, coverages, min_expected_mutation_frequency);
}

void merge(std::vector<GenomicRegion>&& src, std::vector<GenomicRegion>& dst)
{
    const auto itr = utils::append(std::move(src), dst);
    std::inplace_merge(std::begin(dst), itr, std::end(dst));
}

std::vector<GenomicRegion>
get_interesting_hotspots(const GenomicRegion& region,
                         const std::unordered_map<SampleName, CoverageTracker<GenomicRegion>>& interesting_read_tracker,
                         const std::unordered_map<SampleName, CoverageTracker<GenomicRegion>>& tracker,
                         const double min_expected_mutation_frequency)
{
    if (interesting_read_tracker.empty()) {
        return {};
    } else if (interesting_read_tracker.size() == 1) {
        assert(!tracker.empty());
        const auto itr = std::cbegin(interesting_read_tracker);
        assert(tracker.count(itr->first) == 1);
        return get_interesting_hotspots(region, itr->second, tracker.at(itr->first), min_expected_mutation_frequency);
    } else {
        const auto n = size(region);
        std::vector<unsigned> pooled_interesting_coverage(n), pooled_coverage(n);
        std::vector<unsigned> best_sample_interesting_coverage(n), best_sample_coverage(n);
        for (const auto& p : interesting_read_tracker) {
            assert(tracker.count(p.first) == 1);
            const auto sample_coverage = tracker.at(p.first).get(region);
            const auto sample_interesting_coverage = p.second.get(region);
            assert(sample_coverage.size() == n);
            assert(sample_interesting_coverage.size() == n);
            for (std::size_t i {0}; i < sample_coverage.size(); ++i) {
                if (sample_interesting_coverage[i] > 0) {
                    pooled_interesting_coverage[i] += sample_interesting_coverage[i];
                    pooled_coverage[i] += sample_coverage[i];
                    const auto coverage_diff = sample_coverage[i] - sample_interesting_coverage[i];
                    if (coverage_diff > (best_sample_coverage[i] - best_sample_interesting_coverage[i])) {
                        best_sample_interesting_coverage[i] = sample_interesting_coverage[i];
                        best_sample_coverage[i] = sample_coverage[i];
                    }
                }
            }
        }
        auto pooled_hotspots = get_interesting_hotspots(region, pooled_interesting_coverage, pooled_coverage, min_expected_mutation_frequency);
        auto best_sample_hotspots = get_interesting_hotspots(region, best_sample_interesting_coverage, best_sample_coverage, min_expected_mutation_frequency);
        if (!best_sample_hotspots.empty() && best_sample_hotspots != pooled_hotspots) {
            merge(std::move(best_sample_hotspots), pooled_hotspots);
            return extract_covered_regions(pooled_hotspots);
        } else {
            return pooled_hotspots;
        }
    }
}

std::vector<GenomicRegion> AssemblerActiveRegionGenerator::generate(const GenomicRegion& region) const
{
    auto interesting_regions = get_interesting_hotspots(region, interesting_read_coverages_, coverage_tracker_, min_expected_mutation_frequency_);
    if (structual_interesting_) {
        for (const auto& p : clipped_coverage_tracker_) {
            auto deletion_regions = get_deletion_hotspots(region, p.second, read_profile_);
            merge(std::move(deletion_regions), interesting_regions);
        }
    }
    return join(extract_covered_regions(interesting_regions), 30);
}

void AssemblerActiveRegionGenerator::clear() noexcept
{
    free_memory(coverage_tracker_);
    free_memory(interesting_read_coverages_);
    free_memory(clipped_coverage_tracker_);
}

// private methods

namespace {

using NucleotideSequenceIterator = AlignedRead::NucleotideSequence::const_iterator;
using BaseQualityVectorIterator = AlignedRead::BaseQualityVector::const_iterator;

bool is_good_clip(const NucleotideSequenceIterator first_base, const NucleotideSequenceIterator last_base,
                  const BaseQualityVectorIterator first_quality, const AlignedRead::BaseQuality good_base_trigger,
                  const std::size_t min_good_bases)
{
    using boost::make_zip_iterator;
    using Tuple = boost::tuple<char, AlignedRead::BaseQuality>;
    const auto last_quality = std::next(first_quality, std::distance(first_base, last_base));
    return static_cast<std::size_t>(std::count_if(make_zip_iterator(boost::make_tuple(first_base, first_quality)),
                                                  make_zip_iterator(boost::make_tuple(last_base, last_quality)),
                                                  [=] (const Tuple& t) {
                                                      return t.get<0>() != 'N' && t.get<1>() >= good_base_trigger;
                                                  })) >= min_good_bases;
}

} // namespace

bool AssemblerActiveRegionGenerator::is_interesting(const AlignedRead& read) const
{
    using std::cbegin; using std::next;
    using Flag = CigarOperation::Flag;
    const auto& read_sequence = read.sequence();
    auto sequence_itr = cbegin(read_sequence);
    auto base_quality_itr = cbegin(read.base_qualities());
    auto ref_index = mapped_begin(read);
    std::size_t read_index {0};
	std::vector<std::size_t> error_indices {};
    for (const auto& cigar_operation : read.cigar()) {
        const auto op_size = cigar_operation.size();
        switch (cigar_operation.flag()) {
            case Flag::alignmentMatch:
            {
                if (snvs_interesting_ || clustered_interesting_) {
                    const GenomicRegion region {contig_name(read), ref_index, ref_index + op_size};
                    const auto ref_segment = reference_.get().fetch_sequence(region);
					for (std::size_t base_index {0}; base_index < op_size; ++base_index) {
						const auto ref_base = ref_segment[base_index];
						const auto read_base = read.sequence()[read_index + base_index];
						const auto base_quality = read.base_qualities()[read_index + base_index];
						if (ref_base != read_base && ref_base != 'N' && read_base != 'N' && base_quality >= trigger_quality_) {
							if (snvs_interesting_) {
								return true;
							} else { // clustered_interesting_
								error_indices.push_back(read_index + base_index);
							}
						}
					}
                }
                read_index += op_size;
                ref_index += op_size;
                break;
            }
            case Flag::sequenceMatch:
                read_index += op_size;
                ref_index += op_size;
                break;
            case Flag::substitution:
            {
                if (snvs_interesting_ || clustered_interesting_) {
                	GenomicRegion {contig_name(read), ref_index, ref_index + op_size};
					for (std::size_t base_index {0}; base_index < op_size; ++base_index) {
						if (read.base_qualities()[read_index + base_index] >= trigger_quality_) {
							if (snvs_interesting_) {
								return true;
							} else { // clustered_interesting_
								error_indices.push_back(read_index + base_index);
							}
						}
					}
                }
                read_index += op_size;
                ref_index += op_size;
                break;
            }
            case Flag::insertion: return true;
            case Flag::deletion: return true;
            case Flag::softClipped:
            {
                if (indels_interesting_ || structual_interesting_) {
                    if (is_good_clip(next(sequence_itr, read_index), next(sequence_itr, read_index + op_size),
                                     next(base_quality_itr, read_index), trigger_quality_, trigger_clip_size_)) {
                        return true;
                    }
                }
                read_index += op_size;
                ref_index += op_size;
                break;
            }
            case Flag::hardClipped: break;
            case Flag::padding:
                ref_index += op_size;
                break;
            case Flag::skipped:
                ref_index += op_size;
                break;
        }
    }
	if (error_indices.empty()) return false;
	std::vector<unsigned> error_gaps(error_indices.size());
	std::adjacent_difference(std::cbegin(error_indices), std::cend(error_indices), std::begin(error_gaps));
	const auto is_small_gap = [] (auto gap) { return gap <= 10; };
    return std::count_if(std::next(std::cbegin(error_gaps)), std::cend(error_gaps), is_small_gap) >= 3;
}
    
} // namespace coretools
} // namespace octopus
