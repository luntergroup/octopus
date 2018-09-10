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

static constexpr std::array<std::array<double, 31>, 11> enrichment_model
{{
 {1.0, 1.0, 1.05, 1.2, 2.59, 5.34, 9.72, 25.71, 83.84, 136.51, 219.98, 377.14, 611.93, 790.82, 969.93, 1105.94, 1204.84, 1467.31, 1680.1, 1768.04, 1872.59, 2006.16, 2126.59, 2210.81, 2300.07, 2383.73, 2492.18, 2582.96, 2667.07, 2749.69, 2855.54},
 {1.0, 1.0, 1.05, 1.2, 2.59, 5.34, 9.72, 25.71, 83.84, 136.51, 219.98, 377.14, 611.93, 790.82, 969.93, 1105.94, 1204.84, 1467.31, 1680.1, 1768.04, 1872.59, 2006.16, 2126.59, 2210.81, 2300.07, 2383.73, 2492.18, 2582.96, 2667.07, 2749.69, 2855.54},
 {1.0, 1.0, 1.04, 1.81, 5.69, 25.25, 70.14, 143.31, 204.96, 281.79, 344.17, 403.75, 455.72, 502.57, 548.57, 587.21, 630.57, 674.85, 710.06, 735.21, 757.59, 777.2, 795.2, 812.3, 828.7, 844.55, 860.39, 876.56, 892.45, 908.01, 923.77},
 {1.0, 1.0, 1.07, 2.12, 10.45, 36.3, 64.76, 103.63, 141.6, 173.45, 201.76, 224.95, 246.8, 260.35, 271.52, 279.63, 285.49, 290.76, 295.42, 300.08, 304.01, 308.15, 312.51, 317.08, 321.86, 326.86, 332.06, 337.48, 343.11, 348.96, 355.02},
 {1.0, 1.0, 1.16, 5.76, 29.97, 68.27, 109.44, 145.92, 174.82, 198.8, 214.11, 223.12, 229.23, 233.8, 237.5, 240.2, 242.78, 245.54, 248.48, 251.6, 254.9, 258.37, 262.02, 265.86, 269.86, 269.86, 269.86, 269.86, 269.86, 269.86, 269.86},
 {1.0, 1.0, 1.19, 9.19, 34.35, 70.63, 103.19, 146.64, 199.38, 235.41, 244.46, 251.13, 259.76, 262.04, 264.52, 267.19, 270.07, 273.14, 276.41, 279.88, 279.88, 279.88, 279.88, 279.88, 279.88, 279.88, 279.88, 279.88, 279.88, 279.88, 279.88},
 {1.0, 1.0, 1.62, 18.69, 51.78, 80.78, 100.87, 113.06, 118.52, 121.21, 122.88, 123.78, 124.78, 125.88, 127.07, 128.35, 129.73, 129.73, 129.73, 129.73, 129.73, 129.73, 129.73, 129.73, 129.73, 129.73, 129.73, 129.73, 129.73, 129.73, 129.73},
 {1.0, 1.0, 2.9, 23.72, 59.74, 81.39, 89.56, 94.38, 95.83, 96.39, 97.03, 97.74, 98.53, 99.4, 100.34, 100.34, 100.34, 100.34, 100.34, 100.34, 100.34, 100.34, 100.34, 100.34, 100.34, 100.34, 100.34, 100.34, 100.34, 100.34, 100.34},
 {1.0, 1.0, 5.04, 40.4, 73.2, 93.09, 99.93, 102.32, 103.12, 103.73, 104.41, 105.18, 106.03, 106.03, 106.03, 106.03, 106.03, 106.03, 106.03, 106.03, 106.03, 106.03, 106.03, 106.03, 106.03, 106.03, 106.03, 106.03, 106.03, 106.03, 106.03},
 {1.0, 1.0, 6.07, 35.55, 54.35, 60.28, 62.42, 62.69, 63.01, 63.38, 63.8, 64.27, 64.27, 64.27, 64.27, 64.27, 64.27, 64.27, 64.27, 64.27, 64.27, 64.27, 64.27, 64.27, 64.27, 64.27, 64.27, 64.27, 64.27, 64.27, 64.27},
 {1.0, 1.0, 12.83, 52.15, 70.28, 75.32, 76.86, 77.2, 77.6, 78.05, 78.05, 78.05, 78.05, 78.05, 78.05, 78.05, 78.05, 78.05, 78.05, 78.05, 78.05, 78.05, 78.05, 78.05, 78.05, 78.05, 78.05, 78.05, 78.05, 78.05, 78.05}
}};

double calculate_gap_open_rate(const double base_rate, unsigned period, unsigned num_periods) noexcept
{
    if (period == 0 || num_periods == 0) return base_rate;
    static constexpr auto max_period = static_cast<unsigned>(enrichment_model.size());
    period = std::min(period, max_period);
    num_periods = std::min(num_periods, static_cast<unsigned>(enrichment_model[period - 1].size()));
    return std::min(base_rate * enrichment_model[period - 1][num_periods - 1], 1.0);
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
