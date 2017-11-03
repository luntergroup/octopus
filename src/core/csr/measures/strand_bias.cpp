// Copyright (c) 2016 Daniel Cooke
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
#include "basics/aligned_read.hpp"
#include "utils/maths.hpp"
#include "utils/beta_distribution.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

std::unique_ptr<Measure> StrandBias::do_clone() const
{
    return std::make_unique<StrandBias>(*this);
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
DirectionCounts count_directions(const Container& reads)
{
    auto n_fwd = static_cast<unsigned>(std::count_if(std::cbegin(reads), std::cend(reads),
                                                     [] (const auto& read) { return is_forward(read); }));
    return {n_fwd, static_cast<unsigned>(reads.size()) - n_fwd};
}
using DirectionCountVector = std::vector<DirectionCounts>;

auto get_direction_counts(const HaplotypeSupportMap& support, const unsigned prior = 1)
{
    DirectionCountVector result {};
    result.reserve(support.size());
    for (const auto& p : support) {
        result.push_back(count_directions(p.second));
        result.back().forward += prior;
        result.back().reverse += prior;
    }
    return result;
}

double kl_divergence_beta(double a1, double b1, double a2, double b2) noexcept
{
    using boost::math::beta;
    auto a = beta(a2, b2, 1.0);
    auto b = beta(a1, b1, 1.0);
    using boost::math::digamma;
    auto c = (a1 - a2) * digamma(a1);
    auto d = (b1 - b2) * digamma(b1);
    auto e = (a2 - a1 + b2 - b1) * digamma(a1 + b1);
    return a - b + c + d + e;
}

double kl_divergence(const DirectionCounts& lhs, const DirectionCounts& rhs) noexcept
{
    return kl_divergence_beta(lhs.forward, lhs.reverse, rhs.forward, rhs.reverse);
}

double max_pairwise_kl_divergence(const DirectionCountVector& direction_counts) noexcept
{
    double result {0};
    const auto num_counts = direction_counts.size();
    if (num_counts < 2) return result;
    for (std::size_t i {0}; i < num_counts - 1; ++i) {
        for (auto j = i + 1; j < num_counts; ++j) {
            result = std::max(result, kl_divergence(direction_counts[i], direction_counts[j]));
        }
    }
    return result;
}

auto sample_beta(const DirectionCounts& counts, const std::size_t n)
{
    static std::default_random_engine generator {};
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

Measure::ResultType StrandBias::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto assignments = boost::get<ReadAssignments::ResultType>(facets.at("ReadAssignments").get());
    boost::optional<double> result {0};
    for (const auto& p : assignments) {
        if (call.is_heterozygous(p.first)) {
            const auto direction_counts = get_direction_counts(p.second);
            const auto divergence = calculate_max_prob_different(direction_counts, 10'000, 0.25);
            if (result) {
                result = std::max(*result, divergence);
            } else {
                result = divergence;
            }
        }
    }
    return result;
}

std::string StrandBias::do_name() const
{
    return "SB";
}

std::vector<std::string> StrandBias::do_requirements() const
{
    return {"ReadAssignments"};
}

} // namespace csr
} // namespace octopus