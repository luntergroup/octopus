// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "denovo_model.hpp"

#include <memory>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstdint>
#include <cassert>

#include "tandem/tandem.hpp"
#include "basics/phred.hpp"
#include "utils/maths.hpp"
#include "core/types/variant.hpp"

namespace octopus {

namespace {

auto make_flat_hmm_model(const double snv_mutation_rate, const double indel_mutation_rate) noexcept
{
    auto mutation = static_cast<std::int8_t>(probability_to_phred(snv_mutation_rate).score());
    auto gap_open = mutation;
    auto gap_extension = static_cast<std::int8_t>(probability_to_phred(std::min(100 * indel_mutation_rate, 1.0)).score());
    return hmm::FlatGapMutationModel {mutation, gap_open, gap_extension};
}

auto make_gap_open_model(const double indel_mutation_rate, const std::size_t max_repeat_length)
{
    using Penalty = hmm::VariableGapOpenMutationModel::Penalty;
    std::vector<Penalty> result(max_repeat_length);
    for (unsigned i {0}; i < max_repeat_length; ++i) {
        result[i] = static_cast<Penalty>(probability_to_phred(std::min(std::pow(2, i) * indel_mutation_rate, 1.0)).score());
    }
    return result;
}

auto make_gap_extend_model(const double indel_mutation_rate, const std::size_t max_repeat_length)
{
    using Penalty = hmm::VariableGapOpenMutationModel::Penalty;
    std::vector<Penalty> result(max_repeat_length);
    for (unsigned i {0}; i < max_repeat_length; ++i) {
        result[i] = static_cast<Penalty>(probability_to_phred(std::min(std::pow(10, i) * indel_mutation_rate, 1.0)).score());
    }
    return result;
}

} // namespace

DeNovoModel::DeNovoModel(Parameters parameters, std::size_t num_haplotypes_hint, CachingStrategy caching)
: flat_mutation_model_ {make_flat_hmm_model(parameters.snv_mutation_rate, parameters.indel_mutation_rate)}
, repeat_length_gap_open_model_ {make_gap_open_model(parameters.indel_mutation_rate, hmm::min_flank_pad())}
, repeat_length_gap_extend_model_ {make_gap_extend_model(parameters.indel_mutation_rate, hmm::min_flank_pad())}
, min_ln_probability_ {}
, num_haplotypes_hint_ {num_haplotypes_hint}
, haplotypes_ {}
, caching_ {caching}
, gap_open_penalties_ {}
, value_cache_ {}
, address_cache_ {}
, guarded_index_cache_ {}
, unguarded_index_cache_ {}
, padded_given_ {}
, use_unguarded_ {false}
{
    gap_open_penalties_.reserve(1000);
    min_ln_probability_ = 8 * std::log(parameters.snv_mutation_rate);
    if (caching_ == CachingStrategy::address) {
        address_cache_.reserve(num_haplotypes_hint_ * num_haplotypes_hint_);
    } else if (caching == CachingStrategy::value) {
        value_cache_.reserve(num_haplotypes_hint_);
    }
    padded_given_.reserve(1000);
}

void DeNovoModel::prime(std::vector<Haplotype> haplotypes)
{
    constexpr std::size_t max_unguardered {50};
    if (haplotypes.size() <= max_unguardered) {
        unguarded_index_cache_.assign(haplotypes.size(), std::vector<double>(haplotypes.size(), 0));
        for (unsigned target {0}; target < haplotypes.size(); ++target) {
            for (unsigned given {0}; given < haplotypes.size(); ++given) {
                if (target != given) {
                    unguarded_index_cache_[target][given] = evaluate_uncached(haplotypes[target], haplotypes[given]);
                }
            }
        }
        use_unguarded_ = true;
    } else {
        haplotypes_ = std::move(haplotypes);
        guarded_index_cache_.assign(haplotypes_.size(), std::vector<boost::optional<double>>(haplotypes_.size()));
    }
}

void DeNovoModel::unprime() noexcept
{
    haplotypes_.clear();
    haplotypes_.shrink_to_fit();
    guarded_index_cache_.clear();
    guarded_index_cache_.shrink_to_fit();
    unguarded_index_cache_.clear();
    unguarded_index_cache_.shrink_to_fit();
    use_unguarded_ = false;
}

bool DeNovoModel::is_primed() const noexcept
{
    return !guarded_index_cache_.empty() || !unguarded_index_cache_.empty();
}

double DeNovoModel::evaluate(const Haplotype& target, const Haplotype& given) const
{
    switch (caching_) {
        case CachingStrategy::address: return evaluate_address_cache(target, given);
        case CachingStrategy::value: return evaluate_basic_cache(target, given);
        case CachingStrategy::none: return evaluate_uncached(target, given);
        default: return evaluate_uncached(target, given); // to prevent compiler warning
    }
}

double DeNovoModel::evaluate(const unsigned target, const unsigned given) const noexcept
{
    if (use_unguarded_) {
        return unguarded_index_cache_[target][given];
    } else {
        auto& result = guarded_index_cache_[target][given];
        if (!result) {
            if (target != given) {
                result = evaluate_uncached(haplotypes_[target], haplotypes_[given]);
            } else {
                result = 0;
            }
        }
        return *result;
    }
}

// private methods

namespace {

auto get_short_tandem_repeats(const Haplotype::NucleotideSequence& given)
{
    return tandem::extract_exact_tandem_repeats(given, 1, 3);
}

void pad_given(const Haplotype::NucleotideSequence& target, const Haplotype::NucleotideSequence& given,
               std::string& result)
{
    const auto required_size = std::max(target.size(), given.size()) + 2 * hmm::min_flank_pad();
    result.resize(required_size);
    auto itr = std::fill_n(std::begin(result), hmm::min_flank_pad(), 'N');
    itr = std::copy(std::cbegin(given), std::cend(given), itr);
    std::fill(itr, std::end(result), 'N');
}

void pad_given(const Haplotype& target, const Haplotype& given, std::string& result)
{
    pad_given(target.sequence(), given.sequence(), result);
}

auto sequence_length_distance(const Haplotype& lhs, const Haplotype& rhs) noexcept
{
    const auto p = std::minmax({sequence_size(lhs), sequence_size(rhs)});
    return p.second - p.first;
}

bool can_align_with_hmm(const Haplotype& target, const Haplotype& given) noexcept
{
    return sequence_length_distance(target, given) <= hmm::min_flank_pad();
}

template <typename Container>
void rotate(Container& c, const std::size_t n)
{
    std::rotate(std::begin(c), std::next(std::begin(c), n), std::end(c));
}

double hmm_align(const Haplotype::NucleotideSequence& target, const Haplotype::NucleotideSequence& given,
                 const hmm::FlatGapMutationModel& model, const boost::optional<double> min_ln_probability) noexcept
{
    const auto p = hmm::evaluate(target, given, model);
    return min_ln_probability ? std::max(p, *min_ln_probability) : p;
}

double hmm_align(const Haplotype::NucleotideSequence& target, const Haplotype::NucleotideSequence& given,
                 const hmm::VariableGapOpenMutationModel& model, const boost::optional<double> min_ln_probability) noexcept
{
    const auto p = hmm::evaluate(target, given, model);
    return min_ln_probability ? std::max(p, *min_ln_probability) : p;
}

double approx_align(const Haplotype& target, const Haplotype& given, const hmm::FlatGapMutationModel& model,
                    const boost::optional<double> min_ln_probability)
{
    using maths::constants::ln10Div10;
    const auto indel_size = sequence_length_distance(target, given);
    double score = {model.gap_open + model.gap_extend * (static_cast<double>(indel_size) - 1)};
    if (min_ln_probability) {
        const auto max_score = -*min_ln_probability / ln10Div10<>;
        if (score >= max_score) {
            return *min_ln_probability;
        }
    }
    const auto variants = target.difference(given);
    const auto mismatch_size = std::accumulate(std::cbegin(variants), std::cend(variants), 0,
                                               [](auto curr, const auto& variant) noexcept {
                                                   return curr + (!is_indel(variant) ? region_size(variant) : 0);
                                               });
    score += mismatch_size * model.mutation;
    return min_ln_probability ? std::max(-ln10Div10<> * score, *min_ln_probability) : -ln10Div10<> * score;
}

} // namespace

boost::optional<unsigned> DeNovoModel::set_gap_open_penalties(const Haplotype::NucleotideSequence& given) const
{
    const auto repeats = get_short_tandem_repeats(given);
    if (!repeats.empty()) {
        gap_open_penalties_.assign(given.size(), flat_mutation_model_.gap_open);
        const auto max_length = static_cast<unsigned>(repeat_length_gap_open_model_.size());
        unsigned result {0};
        for (const auto& repeat : repeats) {
            const auto num_periods = repeat.period * (repeat.length / repeat.period);
            const auto penalty = repeat_length_gap_open_model_[std::min(num_periods, max_length)];
            std::fill_n(std::next(std::begin(gap_open_penalties_), repeat.pos), repeat.length, penalty);
            result = std::max(repeat.length, result);
        }
        return result;
    } else {
        return boost::none;
    }
}

hmm::VariableGapOpenMutationModel DeNovoModel::make_variable_hmm_model(const unsigned max_repeat_length) const
{
    return {flat_mutation_model_.mutation, gap_open_penalties_, repeat_length_gap_extend_model_[max_repeat_length]};
}

double DeNovoModel::evaluate_uncached(const Haplotype& target, const Haplotype& given) const
{
    if (sequence_size(target) == sequence_size(given)) {
        pad_given(target, given, padded_given_);
        return hmm_align(target.sequence(), padded_given_, flat_mutation_model_, min_ln_probability_);
    } else {
        const auto max_repeat_length = set_gap_open_penalties(given.sequence());
        if (max_repeat_length) {
            if (can_align_with_hmm(target, given)) {
                pad_given(target, given, padded_given_);
                gap_open_penalties_.resize(padded_given_.size(), flat_mutation_model_.gap_open);
                rotate(gap_open_penalties_, hmm::min_flank_pad());
                const auto model = make_variable_hmm_model(*max_repeat_length);
                return hmm_align(target.sequence(), padded_given_, model, min_ln_probability_);
            } else {
                return approx_align(target, given, flat_mutation_model_, min_ln_probability_);
            }
        } else {
            if (can_align_with_hmm(target, given)) {
                pad_given(target, given, padded_given_);
                return hmm_align(target.sequence(), padded_given_, flat_mutation_model_, min_ln_probability_);
            } else {
                return approx_align(target, given, flat_mutation_model_, min_ln_probability_);
            }
        }
    }
}

double DeNovoModel::evaluate_basic_cache(const Haplotype& target, const Haplotype& given) const
{
    const auto target_iter = value_cache_.find(target);
    if (target_iter != std::cend(value_cache_)) {
        const auto given_iter = target_iter->second.find(given);
        if (given_iter != std::cend(target_iter->second)) {
            return given_iter->second;
        }
    }
    const auto result = evaluate_uncached(target, given);
    if (target_iter != std::cend(value_cache_)) {
        target_iter->second.emplace(given, result);
    } else {
        auto p = value_cache_.emplace(std::piecewise_construct,
                                      std::forward_as_tuple(target),
                                      std::forward_as_tuple(num_haplotypes_hint_));
        assert(p.second);
        p.first->second.emplace(given, result);
    }
    return result;
}

double DeNovoModel::evaluate_address_cache(const Haplotype& target, const Haplotype& given) const
{
    const auto p = std::make_pair(std::addressof(target), std::addressof(given));
    const auto itr = address_cache_.find(p);
    if (itr == std::cend(address_cache_)) {
        const auto result = evaluate_uncached(target, given);
        address_cache_.emplace(std::piecewise_construct,
                               std::forward_as_tuple(p),
                               std::forward_as_tuple(result));
        return result;
    } else {
        return itr->second;
    }
}

} // namespace octopus
