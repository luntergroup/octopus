// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "indel_mutation_model.hpp"

#include <utility>
#include <cassert>

#include "utils/maths.hpp"
#include "utils/repeat_finder.hpp"
#include "indel_mutation_model_config.hpp"

namespace octopus {

namespace {

double calculate_gap_open_rate(const double base_rate, unsigned period, unsigned num_periods) noexcept
{
    if (period == 0 || num_periods == 0) return base_rate;
    static constexpr auto max_period = static_cast<unsigned>(enrichment_model.size() - 1);
    static constexpr auto max_periods = static_cast<unsigned>(enrichment_model.front().size() - 1);
    period = std::min(period, max_period);
    num_periods = std::min(num_periods, max_periods);
    return std::min(base_rate * enrichment_model[period][num_periods], 1.0);
}

double calculate_gap_extend_rate(unsigned period, unsigned num_periods, unsigned current_gap, const double gap_open_rate) noexcept
{
    static constexpr auto max_period = static_cast<unsigned>(extension_model.size() - 1);
    static constexpr auto max_periods = static_cast<unsigned>(extension_model.front().size() - 1);
    static constexpr auto max_gap = static_cast<unsigned>(extension_model.front().front().size() - 1);
    period = std::min(period, max_period);
    num_periods = std::min(num_periods, max_periods);
    current_gap = std::min(current_gap, max_gap);
    return std::max(gap_open_rate, extension_model[period][num_periods][current_gap]);
}

} // namespace

IndelMutationModel::IndelMutationModel(Parameters params)
: params_ {std::move(params)}
, indel_repeat_model_ {params_.max_period + 1, std::vector<ModelCell>(params_.max_periodicity + 1)}
{
    for (unsigned period {0}; period <= params_.max_period; ++period) {
        for (unsigned periods {0}; periods <= params_.max_periodicity; ++periods) {
            const auto open_rate = calculate_gap_open_rate(params.indel_mutation_rate, period, periods);
            indel_repeat_model_[period][periods].open = std::min(open_rate, params_.max_open_probability);
            indel_repeat_model_[period][periods].extend.resize(params_.max_indel_length);
            for (unsigned gap {0}; gap < params_.max_indel_length; ++gap) {
                const auto extend_rate = calculate_gap_extend_rate(period, periods, gap, open_rate);
                indel_repeat_model_[period][periods].extend[gap] = std::min(extend_rate, params_.max_extend_probability);
            }
        }
    }
}

namespace {

auto find_short_tandem_repeats(const Haplotype& haplotype)
{
    constexpr unsigned max_repeat_period {5};
    return find_exact_tandem_repeats(haplotype.sequence(), haplotype.mapped_region(), 1, max_repeat_period);
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
        assert(repeat.period() > 0 && repeat.period() <= params_.max_period);
        const auto repeat_offset = static_cast<std::size_t>(begin_distance(haplotype, repeat));
        const auto repeat_len = region_size(repeat);
        const auto num_repeats = static_cast<unsigned>(repeat_len / repeat.period());
        assert(num_repeats > 0);
        const auto& repeat_state = indel_repeat_model_[repeat.period()][std::min(num_repeats, params_.max_periodicity)];
        assert(repeat_offset + repeat_len <= result.gap_open.size());
        assert(repeat_offset + repeat_len <= result.gap_extend.size());
        for (auto pos = repeat_offset; pos < (repeat_offset + repeat_len); ++pos) {
            if (result.gap_open[pos] < repeat_state.open) {
                result.gap_open[pos] = repeat_state.open;
                result.gap_extend[pos] = repeat_state.extend;
            }
        }
    }
    return result;
}

// non-member methods

IndelMutationModel::ContextIndelModel make_indel_model(const Haplotype& context, IndelMutationModel::Parameters params)
{
    IndelMutationModel model {params};
    return model.evaluate(context);
}

IndelMutationModel::Probability
calculate_indel_probability(const IndelMutationModel::ContextIndelModel& model, const std::size_t pos, const std::size_t length) noexcept
{
    assert(length > 0 && pos < model.gap_open.size() && pos < model.gap_extend.size());
    const auto& gap_extend_state = model.gap_extend[pos];
    const auto max_gap_length = gap_extend_state.size() - 1;
    const auto gap_extend_itr = std::next(std::cbegin(gap_extend_state));
    return std::accumulate(gap_extend_itr, std::next(gap_extend_itr, std::min(length - 1, max_gap_length)),
                           model.gap_open[pos], std::multiplies<> {});
}

} // namespace octopus
