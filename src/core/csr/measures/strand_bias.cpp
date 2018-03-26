// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "strand_bias.hpp"

#include <algorithm>
#include <iterator>
#include <random>
#include <functional>
#include <cmath>
#include <cassert>

#include <boost/variant.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>

#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "basics/aligned_read.hpp"
#include "utils/maths.hpp"
#include "utils/beta_distribution.hpp"
#include "utils/genotype_reader.hpp"
#include "../facets/samples.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

StrandBias::StrandBias(const double critical_value)
: min_medium_trigger_ {critical_value / 2}
, min_big_trigger_ {critical_value / 8}
, critical_resample_lb_ {0.995 * critical_value}
, critical_resample_ub_ {1.005 * critical_value}
, use_resampling_ {true}
{}

std::unique_ptr<Measure> StrandBias::do_clone() const
{
    return std::make_unique<StrandBias>(*this);
}

namespace {

bool is_canonical(const VcfRecord::NucleotideSequence& allele) noexcept
{
    const static VcfRecord::NucleotideSequence deleted_allele {vcfspec::deletedBase};
    return !(allele == vcfspec::missingValue || allele == deleted_allele);
}

bool has_called_alt_allele(const VcfRecord& call, const VcfRecord::SampleName& sample)
{
    if (!call.has_genotypes()) return true;
    const auto& genotype = get_genotype(call, sample);
    return std::any_of(std::cbegin(genotype), std::cend(genotype),
                       [&] (const auto& allele) { return allele != call.ref() && is_canonical(allele); });
}

bool is_evaluable(const VcfRecord& call, const VcfRecord::SampleName& sample)
{
    return has_called_alt_allele(call, sample) && call.is_heterozygous(sample);
}

bool is_forward(const AlignedRead& read) noexcept
{
    return read.direction() == AlignedRead::Direction::forward;
}

struct DirectionCounts
{
    unsigned forward, reverse;
};

template <typename Container>
DirectionCounts count_directions(const Container& reads, const GenomicRegion& call_region)
{
    unsigned n_forward {0}, n_reverse {0};
    for (const auto& read : reads) {
        if (overlaps(read.get(), call_region)) {
            if (is_forward(read)) {
                ++n_forward;
            } else {
                ++n_reverse;
            }
        }
    }
    return {n_forward, n_reverse};
}

using DirectionCountVector = std::vector<DirectionCounts>;

auto get_direction_counts(const AlleleSupportMap& support, const GenomicRegion& call_region, const unsigned prior = 1)
{
    DirectionCountVector result {};
    result.reserve(support.size());
    for (const auto& p : support) {
        result.push_back(count_directions(p.second, call_region));
        result.back().forward += prior;
        result.back().reverse += prior;
    }
    return result;
}

auto sample_beta(const DirectionCounts& counts, const std::size_t n)
{
    static thread_local std::mt19937 generator {42};
    std::beta_distribution<> beta {static_cast<double>(counts.forward), static_cast<double>(counts.reverse)};
    std::vector<double> result(n);
    std::generate_n(std::begin(result), n, [&] () { return beta(generator); });
    return result;
}

auto generate_beta_samples(const DirectionCountVector& direction_counts, const std::size_t num_samples)
{
    std::vector<std::vector<double>> result {};
    result.reserve(direction_counts.size());
    for (const auto& counts : direction_counts) {
        result.push_back(sample_beta(counts, num_samples));
    }
    return result;
}

double estimate_prob_different(const std::vector<double>& lhs, const std::vector<double>& rhs,
                               const double min_diff)
{
    assert(lhs.size() == rhs.size());
    std::vector<double> diffs(lhs.size());
    std::transform(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs), std::begin(diffs), std::minus<> {});
    auto n_diffs = std::count_if(std::cbegin(diffs), std::cend(diffs), [=] (auto diff) { return std::abs(diff) > min_diff; });
    return static_cast<double>(n_diffs) / diffs.size();
}

double calculate_max_prob_different(const DirectionCountVector& direction_counts, const std::size_t num_samples,
                                    const double min_diff)
{
    const auto num_counts = direction_counts.size();
    if (num_counts < 2) return 0;
    const auto samples = generate_beta_samples(direction_counts, num_samples);
    double result {0};
    for (std::size_t i {0}; i < num_counts - 1; ++i) {
        for (auto j = i + 1; j < num_counts; ++j) {
            result = std::max(result, estimate_prob_different(samples[i], samples[j], min_diff));
        }
    }
    return result;
}

} // namespace

Measure::ResultType StrandBias::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).support;
    std::vector<boost::optional<double>> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        boost::optional<double> sample_result {};
        if (is_evaluable(call, sample)) {
            std::vector<Allele> alleles; bool has_ref;
            std::tie(alleles, has_ref) = get_called_alleles(call, sample, true);
            assert(!alleles.empty());
            const auto sample_allele_support = compute_allele_support(alleles, assignments.at(sample));
            const auto direction_counts = get_direction_counts(sample_allele_support, mapped_region(call));
            double prob;
            if (use_resampling_) {
                prob = calculate_max_prob_different(direction_counts, small_sample_size_, min_difference_);
                if (prob >= min_big_trigger_) {
                    prob = calculate_max_prob_different(direction_counts, big_sample_size_, min_difference_);
                } else if (prob >= min_medium_trigger_) {
                    prob = calculate_max_prob_different(direction_counts, medium_sample_size_, min_difference_);
                    if (prob >= min_big_trigger_) {
                        prob = calculate_max_prob_different(direction_counts, big_sample_size_, min_difference_);
                    }
                }
                if (prob > critical_resample_lb_ && prob < critical_resample_ub_) {
                    prob = calculate_max_prob_different(direction_counts, very_big_sample_size, min_difference_);
                }
            } else {
                prob = calculate_max_prob_different(direction_counts, big_sample_size_, min_difference_);
            }
            sample_result = prob;
        }
        result.push_back(sample_result);
    }
    return result;
}

Measure::ResultCardinality StrandBias::do_cardinality() const noexcept
{
    return ResultCardinality::num_samples;
}

std::string StrandBias::do_name() const
{
    return "SB";
}

std::string StrandBias::do_describe() const
{
    return "Strand bias of reads based on haplotype support";
}

std::vector<std::string> StrandBias::do_requirements() const
{
    return {"Samples", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus