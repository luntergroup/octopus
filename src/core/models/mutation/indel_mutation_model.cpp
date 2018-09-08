// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "indel_mutation_model.hpp"

#include <utility>
#include <array>
#include <cassert>

#include "utils/maths.hpp"
#include "utils/repeat_finder.hpp"

namespace octopus {

namespace {

// The following enrichment model is derived from Montgomery et al. Genome Research. 2013.

static constexpr std::array<std::array<double, 30>, 10> enrichment_model
{{
{4.6, 4.7, 4.8, 6.2, 8.9, 13.2, 28.9, 85.8, 136.7, 216.7, 366.2, 587.3, 751.3, 912.1, 1029.1, 1109.0, 1334.6, 1509.6, 1568.8, 1639.9, 1733.1, 1811.3, 1855.8, 1902.0, 1941.0, 1997.3, 2036.8, 2068.5, 2096.8, 2140.3},
{4.6, 4.6, 5.4, 9.3, 28.7, 73.0, 144.9, 204.6, 278.4, 337.0, 391.8, 438.3, 478.8, 517.4, 548.1, 582.1, 615.8, 640.1, 654.5, 665.6, 673.6, 679.6, 684.1, 687.6, 690.0, 691.9, 693.6, 694.6, 694.8, 694.8},
{4.6, 4.7, 5.7, 14.0, 39.7, 67.7, 105.8, 142.5, 172.7, 199.0, 219.9, 239.0, 249.8, 257.9, 262.9, 265.5, 267.4, 268.4, 269.3, 269.3, 269.3, 269.3, 269.3, 269.3, 269.3, 269.3, 269.3, 269.3, 269.3, 269.3},
{4.6, 4.8, 9.4, 33.5, 71.4, 112.0, 147.4, 175.1, 197.4, 211.0, 218.1, 222.2, 224.7, 226.1, 226.3, 226.3, 226.3, 226.3, 226.3, 226.3, 226.3, 226.3, 226.3, 226.3, 226.3, 226.3, 226.3, 226.3, 226.3, 226.3},
{4.6, 4.8, 12.8, 37.8, 73.8, 105.8, 148.2, 199.1, 233.1, 240.4, 245.1, 251.4, 251.4, 251.4, 251.4, 251.4, 251.4, 251.4, 251.4, 251.4, 251.4, 251.4, 251.4, 251.4, 251.4, 251.4, 251.4, 251.4, 251.4, 251.4},
{4.6, 5.2, 22.3, 55.2, 83.9, 103.5, 115.1, 119.8, 121.8, 122.6, 122.6, 122.6, 122.6, 122.6, 122.6, 122.6, 122.6, 122.6, 122.6, 122.6, 122.6, 122.6, 122.6, 122.6, 122.6, 122.6, 122.6, 122.6, 122.6, 122.6},
{4.6, 6.5, 27.3, 63.1, 84.5, 92.3, 96.6, 97.6, 97.6, 97.6, 97.6, 97.6, 97.6, 97.6, 97.6, 97.6, 97.6, 97.6, 97.6, 97.6, 97.6, 97.6, 97.6, 97.6, 97.6, 97.6, 97.6, 97.6, 97.6, 97.6},
{4.6, 8.6, 43.9, 76.5, 96.1, 102.5, 104.5, 104.7, 104.7, 104.7, 104.7, 104.7, 104.7, 104.7, 104.7, 104.7, 104.7, 104.7, 104.7, 104.7, 104.7, 104.7, 104.7, 104.7, 104.7, 104.7, 104.7, 104.7, 104.7, 104.7},
{4.6, 9.7, 39.1, 57.8, 63.5, 65.4, 65.4, 65.4, 65.4, 65.4, 65.4, 65.4, 65.4, 65.4, 65.4, 65.4, 65.4, 65.4, 65.4, 65.4, 65.4, 65.4, 65.4, 65.4, 65.4, 65.4, 65.4, 65.4, 65.4, 65.4},
{4.6, 16.4, 55.7, 73.6, 78.4, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7}
}};

double calculate_gap_open_rate(const double base_rate, unsigned period, unsigned num_periods) noexcept
{
    if (period == 0 || num_periods == 0) return base_rate;
    static constexpr auto max_period = static_cast<unsigned>(enrichment_model.size());
    period = std::min(period, max_period);
    num_periods = std::min(num_periods, static_cast<unsigned>(enrichment_model[period - 1].size()));
    return std::min(base_rate * enrichment_model[period - 1][num_periods - 1] / enrichment_model[0][0], 1.0);
}

double calculate_gap_extend_rate(const double base_rate, const unsigned period, const unsigned num_periods, const double gap_open_rate)
{
    return std::max(1'000 * gap_open_rate, 0.7);
}

} // namespace

IndelMutationModel::IndelMutationModel(Parameters params)
: params_ {std::move(params)}
, indel_repeat_model_ {params_.max_period + 1, std::vector<ModelCell>(params_.max_periodicity + 1)}
{
    for (unsigned period {0}; period <= params_.max_period; ++period) {
        for (unsigned n {0}; n <= params_.max_periodicity; ++n) {
            const auto open_rate = calculate_gap_open_rate(params.indel_mutation_rate, period, n);
            indel_repeat_model_[period][n].open = std::min(open_rate, params_.max_open_probability);
            const auto extend_rate = calculate_gap_extend_rate(params.indel_mutation_rate, period, n, open_rate);
            indel_repeat_model_[period][n].extend = std::min(extend_rate, params_.max_extend_probability);
        }
    }
}

namespace {

auto find_short_tandem_repeats(const Haplotype& haplotype)
{
    constexpr unsigned max_repeat_period {5};
    return find_exact_tandem_repeats(haplotype.sequence(), haplotype.mapped_region(), 1, max_repeat_period);
}

template <typename FordwardIt, typename Tp>
auto fill_if_greater(FordwardIt first, FordwardIt last, const Tp& value)
{
    return std::transform(first, last, first, [&] (const auto& x) { return std::max(x, value); });
}

template <typename FordwardIt, typename Tp>
auto fill_n_if_greater(FordwardIt first, std::size_t n, const Tp& value)
{
    return fill_if_greater(first, std::next(first, n), value);
}

} // namespace

IndelMutationModel::ContextIndelModel IndelMutationModel::evaluate(const Haplotype& haplotype) const
{
    const auto repeats = find_short_tandem_repeats(haplotype);
    ContextIndelModel result {};
    const auto haplotype_len = sequence_size(haplotype);
    const auto& base_probabilities = indel_repeat_model_[0][0];
    result.gap_open.resize(haplotype_len, base_probabilities.open);
    result.gap_extend.resize(haplotype_len, base_probabilities.extend);
    for (const auto& repeat : repeats) {
        assert(repeat.period > 0 && repeat.period <= params_.max_period);
        const auto repeat_offset = static_cast<std::size_t>(begin_distance(haplotype, repeat));
        const auto repeat_len = region_size(repeat);
        const auto num_repeats = static_cast<unsigned>(repeat_len / repeat.period);
        assert(num_repeats > 0);
        const auto& repeat_state = indel_repeat_model_[repeat.period][std::min(num_repeats, params_.max_periodicity)];
        assert(repeat_offset + repeat_len <= result.gap_open.size());
        fill_n_if_greater(std::next(std::begin(result.gap_open), repeat_offset), repeat_len, repeat_state.open);
        assert(repeat_offset + repeat_len <= result.gap_extend.size());
        fill_n_if_greater(std::next(std::begin(result.gap_extend), repeat_offset), repeat_len, repeat_state.extend);
    }
    return result;
}

// non-member methods

IndelMutationModel::ContextIndelModel make_indel_model(const Haplotype& context, IndelMutationModel::Parameters params)
{
    IndelMutationModel model {params};
    return model.evaluate(context);
}

} // namespace octopus
