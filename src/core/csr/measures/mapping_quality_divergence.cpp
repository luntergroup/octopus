// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "mapping_quality_divergence.hpp"

#include <algorithm>
#include <iterator>
#include <cmath>
#include <limits>
#include <numeric>

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include "basics/aligned_read.hpp"
#include "core/types/allele.hpp"
#include "io/variant/vcf_record.hpp"
#include "utils/maths.hpp"
#include "../facets/samples.hpp"
#include "../facets/alleles.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string MappingQualityDivergence::name_ = "MQD";

std::unique_ptr<Measure> MappingQualityDivergence::do_clone() const
{
    return std::make_unique<MappingQualityDivergence>(*this);
}

Measure::ValueType MappingQualityDivergence::get_value_type() const
{
    return int {};
}

namespace {

using MappingQualityVector = std::vector<AlignedRead::MappingQuality>;

auto extract_mapping_qualities(const ReadRefSupportSet& reads)
{
    MappingQualityVector result(reads.size());
    std::transform(std::cbegin(reads), std::cend(reads), std::begin(result),
                   [](const AlignedRead& read) { return read.mapping_quality(); });
    return result;
}

auto extract_mapping_qualities(const std::vector<Allele>& alleles, const AlleleSupportMap& support, const bool drop_empty = true)
{
    std::vector<MappingQualityVector> result {};
    result.reserve(alleles.size());
    for (const auto& allele : alleles) {
        const auto support_set_itr = support.find(allele);
        if (support_set_itr != std::cend(support) && (!support_set_itr->second.empty() || !drop_empty)) {
            result.push_back(extract_mapping_qualities(support_set_itr->second));
        }
    }
    return result;
}

using MappingQualityHistogram = std::vector<unsigned>;

auto make_histogram(const MappingQualityVector& mapping_qualities, const unsigned prior = 1)
{
    using Q = AlignedRead::MappingQuality;
    static constexpr unsigned max_qualities {std::numeric_limits<Q>::max() - std::numeric_limits<Q>::min()};
    MappingQualityHistogram result(max_qualities, prior);
    for (auto quality : mapping_qualities) {
        ++result[quality];
    }
    return result;
}

auto convert_to_probabilities(const MappingQualityHistogram& counts)
{
    const auto norm = std::accumulate(std::cbegin(counts), std::cend(counts), 0.0);
    std::vector<double> result(counts.size());
    std::transform(std::cbegin(counts), std::cend(counts), std::begin(result),
                   [norm] (const auto c) { return static_cast<double>(c) / norm; });
    return result;
}

auto kl_divergence(const std::vector<double>& p, const std::vector<double>& q) noexcept
{
    return std::inner_product(std::cbegin(p), std::cend(p), std::cbegin(q), 0.0,
                              std::plus<> {}, [] (const auto a, const auto b) {
        return a * std::log(a / b);
    });
}

auto symmetric_kl_divergence(const std::vector<double>& p, const std::vector<double>& q) noexcept
{
    return kl_divergence(p, q) + kl_divergence(q, p);
}

double calculate_max_pairwise_kl_divergence(const std::vector<std::vector<double>>& mq_distributions)
{
    double result {0};
    const auto num_distributions = mq_distributions.size();
    if (num_distributions < 2) return result;
    for (std::size_t i {0}; i < num_distributions - 1; ++i) {
        for (auto j = i + 1; j < num_distributions; ++j) {
            result = std::max(result, symmetric_kl_divergence(mq_distributions[i], mq_distributions[j]));
        }
    }
    return result;
}

double calculate_max_pairwise_kl_divergence(const std::vector<MappingQualityVector>& mapping_qualities)
{
    std::vector<std::vector<double>> mq_distributions {};
    mq_distributions.reserve(mapping_qualities.size());
    for (const auto& mqs : mapping_qualities) {
        mq_distributions.push_back(convert_to_probabilities(make_histogram(mqs)));
    }
    return calculate_max_pairwise_kl_divergence(mq_distributions);
}

auto calculate_median_mapping_qualities(const std::vector<MappingQualityVector>& mapping_qualities)
{
    std::vector<AlignedRead::MappingQuality> result(mapping_qualities.size());
    std::transform(std::cbegin(mapping_qualities), std::cend(mapping_qualities), std::begin(result),
                   [] (const auto& mqs) { return maths::median(mqs); });
    return result;
}

auto calculate_max_pairwise_difference(const std::vector<AlignedRead::MappingQuality>& mapping_qualities)
{
    assert(!mapping_qualities.empty());
    const auto p = std::minmax_element(std::cbegin(mapping_qualities), std::cend(mapping_qualities));
    return *p.second - *p.first;
}

auto max_pairwise_median_mapping_quality_difference(const std::vector<MappingQualityVector>& mapping_qualities)
{
    const auto median_mapping_qualities = calculate_median_mapping_qualities(mapping_qualities);
    return calculate_max_pairwise_difference(median_mapping_qualities);
}

} // namespace

Measure::ResultType MappingQualityDivergence::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& alleles = get_value<Alleles>(facets.at("Alleles"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).alleles;
    Array<Optional<ValueType>> result(samples.size());
    for (std::size_t s {0}; s < samples.size(); ++s) {
        const auto& sample = samples[s];
        if (call.is_heterozygous(sample)) {
            const auto mapping_qualities = extract_mapping_qualities(get_called(alleles, call, sample), assignments.at(sample));
            if (!mapping_qualities.empty()) {
                    result[s] = max_pairwise_median_mapping_quality_difference(mapping_qualities);
            }
        }
    }
    return result;
}

Measure::ResultCardinality MappingQualityDivergence::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& MappingQualityDivergence::do_name() const
{
    return name_;
}

std::string MappingQualityDivergence::do_describe() const
{
    return "Maximum pairwise difference in median mapping qualities of reads supporting each haplotype";
}

std::vector<std::string> MappingQualityDivergence::do_requirements() const
{
    return {"Samples", "Alleles", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
