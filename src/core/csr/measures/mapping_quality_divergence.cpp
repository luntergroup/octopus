// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "mapping_quality_divergence.hpp"

#include <algorithm>
#include <iterator>
#include <cmath>
#include <limits>
#include <numeric>

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "basics/aligned_read.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

std::unique_ptr<Measure> MappingQualityDivergence::do_clone() const
{
    return std::make_unique<MappingQualityDivergence>(*this);
}

namespace {

using MappingQualityVector = std::vector<AlignedRead::MappingQuality>;

auto extract_mapping_qualities(const std::vector<AlignedRead>& reads)
{
    MappingQualityVector result(reads.size());
    std::transform(std::cbegin(reads), std::cend(reads), std::begin(result),
                   [](const auto& read) { return read.mapping_quality(); });
    return result;
}

auto extract_mapping_qualities(const HaplotypeSupportMap& support)
{
    std::vector<MappingQualityVector> result{};
    result.reserve(support.size());
    for (const auto& p : support) {
        result.push_back(extract_mapping_qualities(p.second));
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

} // namespace

Measure::ResultType MappingQualityDivergence::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto assignments = boost::get<ReadAssignments::ResultType>(facets.at("ReadAssignments").get());
    double max_divergence {0};
    for (const auto& p : assignments) {
        if (call.is_heterozygous(p.first)) {
            auto mapping_qualities = extract_mapping_qualities(p.second);
            max_divergence = std::max(max_divergence, calculate_max_pairwise_kl_divergence(mapping_qualities));
        }
    }
    return max_divergence;
}

std::string MappingQualityDivergence::do_name() const
{
    return "MQD";
}

std::vector<std::string> MappingQualityDivergence::do_requirements() const
{
    return {"ReadAssignments"};
}

} // namespace csr
} // namespace octopus
