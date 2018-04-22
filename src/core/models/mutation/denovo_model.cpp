// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "denovo_model.hpp"

#include <memory>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <cassert>

#include "tandem/tandem.hpp"
#include "basics/phred.hpp"
#include "utils/maths.hpp"
#include "core/types/variant.hpp"

namespace octopus {

namespace {

auto make_flat_hmm_model(const double snv_mutation_rate, const double indel_mutation_rate,
                         const double max_rate = 0.5) noexcept
{
    auto mutation_penalty = static_cast<std::int8_t>(probability_to_phred(snv_mutation_rate).score());
    auto gap_open_penalty = mutation_penalty;
    auto gap_extension_penalty = static_cast<std::int8_t>(probability_to_phred(std::min(10'000 * indel_mutation_rate, max_rate)).score());
    return hmm::FlatGapMutationModel {mutation_penalty, gap_open_penalty, gap_extension_penalty};
}

} // namespace

DeNovoModel::DeNovoModel(Parameters parameters, std::size_t num_haplotypes_hint, CachingStrategy caching)
: flat_mutation_model_ {make_flat_hmm_model(parameters.snv_mutation_rate, parameters.indel_mutation_rate)}
, indel_model_ {{parameters.indel_mutation_rate}}
, min_ln_probability_ {}
, num_haplotypes_hint_ {num_haplotypes_hint}
, haplotypes_ {}
, caching_ {caching}
, gap_open_penalties_ {}
, gap_extend_penalties_ {}
, gap_model_index_cache_ {}
, value_cache_ {}
, address_cache_ {}
, guarded_index_cache_ {}
, unguarded_index_cache_ {}
, padded_given_ {}
, use_unguarded_ {false}
{
    gap_open_penalties_.reserve(1000);
    if (caching_ == CachingStrategy::address) {
        address_cache_.reserve(num_haplotypes_hint_ * num_haplotypes_hint_);
    } else if (caching == CachingStrategy::value) {
        value_cache_.reserve(num_haplotypes_hint_);
    }
    padded_given_.reserve(1000);
}

void DeNovoModel::prime(std::vector<Haplotype> haplotypes)
{
    if (is_primed()) throw std::runtime_error {"DeNovoModel: already primed"};
    constexpr std::size_t max_unguardered {50};
    haplotypes_ = std::move(haplotypes);
    gap_model_index_cache_.resize(haplotypes_.size());
    if (haplotypes_.size() <= max_unguardered) {
        unguarded_index_cache_.assign(haplotypes_.size(), std::vector<double>(haplotypes_.size(), 0));
        for (unsigned target {0}; target < haplotypes_.size(); ++target) {
            for (unsigned given {0}; given < haplotypes_.size(); ++given) {
                if (target != given) {
                    unguarded_index_cache_[target][given] = evaluate_uncached(target, given);
                }
            }
        }
        use_unguarded_ = true;
    } else {
        guarded_index_cache_.assign(haplotypes_.size(), std::vector<boost::optional<double>>(haplotypes_.size()));
    }
}

void DeNovoModel::unprime() noexcept
{
    haplotypes_.clear();
    haplotypes_.shrink_to_fit();
    gap_open_penalties_.clear();
    gap_open_penalties_.shrink_to_fit();
    gap_extend_penalties_.clear();
    gap_extend_penalties_.shrink_to_fit();
    gap_model_index_cache_.clear();
    gap_model_index_cache_.shrink_to_fit();
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

double DeNovoModel::evaluate(const unsigned target, const unsigned given) const
{
    if (use_unguarded_) {
        return unguarded_index_cache_[target][given];
    } else {
        auto& result = guarded_index_cache_[target][given];
        if (!result) {
            if (target != given) {
                result = evaluate_uncached(target, given);
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
    constexpr unsigned max_repeat_period {3};
    return tandem::extract_exact_tandem_repeats(given, 1, max_repeat_period);
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
void rotate_right(Container& c, const std::size_t n)
{
    assert(c.size() > n);
    std::rotate(std::rbegin(c), std::next(std::rbegin(c), n), std::rend(c));
}

double calculate_score(const Variant& variant, const Haplotype& context,
                       const hmm::VariableGapExtendMutationModel& mutation_model)
{
    if (is_indel(variant)) {
        assert(contains(context, variant));
        const auto offset = static_cast<std::size_t>(begin_distance(context, variant));
        const auto indel_length = indel_size(variant);
        assert(offset < mutation_model.gap_open.size());
        constexpr decltype(indel_length) max_indel_length {50};
        return mutation_model.gap_open[offset] + (std::min(indel_length, max_indel_length) - 1) * mutation_model.gap_extend[offset];
    } else {
        return mutation_model.mutation;
    }
}

double approx_align(const Haplotype& target, const Haplotype& given,
                    const hmm::VariableGapExtendMutationModel& mutation_model)
{
    double score {0};
    const auto variants = target.difference(given);
    for (const auto& variant : variants) {
        score += calculate_score(variant, given, mutation_model);
    }
    return score * -maths::constants::ln10Div10<>;
}

void set_penalties(const IndelMutationModel::ContextIndelModel::ProbabilityVector& probabilities,
                  hmm::VariableGapExtendMutationModel::PenaltyVector& penalties)
{
    assert(probabilities.size() == penalties.size());
    std::transform(std::cbegin(probabilities), std::cend(probabilities), std::begin(penalties),
                   [] (auto rate) {
                       static constexpr double max_penalty {127};
                       return std::max(std::min(std::log(rate) / -maths::constants::ln10Div10<>, max_penalty), 1.0);
                   });
}

} // namespace

void DeNovoModel::set_gap_penalties(const Haplotype& given) const
{
    const auto contextual_indel_model = indel_model_.evaluate(given);
    const auto num_bases = sequence_size(given);
    assert(contextual_indel_model.gap_open.size() == num_bases);
    gap_open_penalties_.resize(num_bases);
    set_penalties(contextual_indel_model.gap_open, gap_open_penalties_);
    assert(contextual_indel_model.gap_extend.size() == num_bases);
    gap_extend_penalties_.resize(num_bases);
    set_penalties(contextual_indel_model.gap_extend, gap_extend_penalties_);
}

void DeNovoModel::set_gap_penalties(const unsigned given) const
{
    assert(given < gap_model_index_cache_.size());
    auto& cached_result = gap_model_index_cache_[given];
    if (cached_result) {
        gap_open_penalties_   = cached_result->first;
        gap_extend_penalties_ = cached_result->second;
    } else {
        set_gap_penalties(haplotypes_[given]);
        cached_result = std::make_pair(gap_open_penalties_, gap_extend_penalties_);
    }
}

hmm::VariableGapExtendMutationModel DeNovoModel::make_hmm_model_from_cache() const
{
    return {flat_mutation_model_.mutation, gap_open_penalties_, gap_extend_penalties_};
}

double DeNovoModel::evaluate_uncached(const Haplotype& target, const Haplotype& given, const bool gap_penalties_cached) const
{
    if (!gap_penalties_cached) set_gap_penalties(given);
    const auto mutation_model = make_hmm_model_from_cache();
    double result;
    if (can_align_with_hmm(target, given)) {
        pad_given(target, given, padded_given_);
        gap_open_penalties_.resize(padded_given_.size(), mutation_model.mutation);
        rotate_right(gap_open_penalties_, hmm::min_flank_pad());
        gap_extend_penalties_.resize(padded_given_.size(), mutation_model.mutation);
        rotate_right(gap_extend_penalties_, hmm::min_flank_pad());
        result = hmm::evaluate(target.sequence(), padded_given_, mutation_model);
    } else {
        result = approx_align(target, given, mutation_model);
    }
    return min_ln_probability_ ? std::max(result, *min_ln_probability_) : result;
}

double DeNovoModel::evaluate_uncached(const unsigned target_idx, const unsigned given_idx) const
{
    set_gap_penalties(given_idx);
    return evaluate_uncached(haplotypes_[target_idx], haplotypes_[given_idx], true);
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
