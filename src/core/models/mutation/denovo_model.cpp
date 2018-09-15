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

#include "basics/phred.hpp"
#include "utils/maths.hpp"
#include "core/types/variant.hpp"

namespace octopus {

namespace {

std::int8_t probability_to_penalty(const double probability) noexcept
{
    static constexpr double max_penalty {127};
    return std::max(std::min(std::log(probability) / -maths::constants::ln10Div10<>, max_penalty), 1.0);
}

} // namespace

DeNovoModel::DeNovoModel(Parameters parameters, std::size_t num_haplotypes_hint, CachingStrategy caching)
: params_ {parameters}
, snv_penalty_ {probability_to_penalty(params_.snv_mutation_rate)}
, indel_model_ {{params_.indel_mutation_rate}}
, min_ln_probability_ {}
, num_haplotypes_hint_ {num_haplotypes_hint}
, haplotypes_ {}
, caching_ {caching}
, local_indel_model_ {}
, gap_model_index_cache_ {}
, value_cache_ {}
, address_cache_ {}
, guarded_index_cache_ {}
, unguarded_index_cache_ {}
, padded_given_ {}
, use_unguarded_ {false}
{
    if (caching_ == CachingStrategy::address) {
        address_cache_.reserve(num_haplotypes_hint_ * num_haplotypes_hint_);
    } else if (caching == CachingStrategy::value) {
        value_cache_.reserve(num_haplotypes_hint_);
    }
    padded_given_.reserve(1000);
}

DeNovoModel::Parameters DeNovoModel::parameters() const
{
    return params_;
}

void DeNovoModel::prime(std::vector<Haplotype> haplotypes)
{
    if (is_primed()) throw std::runtime_error {"DeNovoModel: already primed"};
    constexpr std::size_t max_unguardered {50};
    haplotypes_ = std::move(haplotypes);
    gap_model_index_cache_.resize(haplotypes_.size());
    if (haplotypes_.size() <= max_unguardered) {
        unguarded_index_cache_.assign(haplotypes_.size(), std::vector<LogProbability>(haplotypes_.size(), 0));
        for (unsigned target {0}; target < haplotypes_.size(); ++target) {
            for (unsigned given {0}; given < haplotypes_.size(); ++given) {
                if (target != given) {
                    unguarded_index_cache_[target][given] = evaluate_uncached(target, given);
                }
            }
        }
        use_unguarded_ = true;
    } else {
        guarded_index_cache_.assign(haplotypes_.size(), std::vector<boost::optional<LogProbability>>(haplotypes_.size()));
    }
}

void DeNovoModel::unprime() noexcept
{
    haplotypes_.clear();
    haplotypes_.shrink_to_fit();
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

DeNovoModel::LogProbability DeNovoModel::evaluate(const Haplotype& target, const Haplotype& given) const
{
    switch (caching_) {
        case CachingStrategy::address: return evaluate_address_cache(target, given);
        case CachingStrategy::value: return evaluate_basic_cache(target, given);
        case CachingStrategy::none: return evaluate_uncached(target, given);
        default: return evaluate_uncached(target, given); // to prevent compiler warning
    }
}

DeNovoModel::LogProbability DeNovoModel::evaluate(const unsigned target, const unsigned given) const
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
    return sequence_length_distance(target, given) < hmm::min_flank_pad();
}

double calculate_log_probability(const Variant& variant, const Haplotype& context,
                                 const double snv_penalty, const IndelMutationModel::ContextIndelModel& indel_model)
{
    if (is_indel(variant)) {
        assert(contains(context, variant));
        const auto offset = static_cast<std::size_t>(begin_distance(context, variant));
        return std::log(calculate_indel_probability(indel_model, offset, indel_size(variant)));
    } else {
        return snv_penalty * -maths::constants::ln10Div10<>;
    }
}

double calculate_approx_log_probability(const Haplotype& target, const Haplotype& given,
                                        const double snv_penalty, const IndelMutationModel::ContextIndelModel& indel_model)
{
    double score {0};
    const auto variants = target.difference(given);
    for (const auto& variant : variants) {
        score += calculate_log_probability(variant, given, snv_penalty, indel_model);
    }
    return score ;
}

void set_penalties(const IndelMutationModel::ContextIndelModel& indel_model,
                   hmm::VariableGapExtendMutationModel::PenaltyVector& open_penalties,
                   hmm::VariableGapExtendMutationModel::PenaltyVector& extend_penalties)
{
    assert(indel_model.gap_open.size() + 2 * hmm::min_flank_pad() == open_penalties.size());
    std::transform(std::cbegin(indel_model.gap_open), std::cend(indel_model.gap_open),
                   std::next(std::begin(open_penalties), hmm::min_flank_pad()),
                   probability_to_penalty);
    assert(indel_model.gap_extend.size() + 2 * hmm::min_flank_pad() == extend_penalties.size());
    std::transform(std::cbegin(indel_model.gap_extend), std::cend(indel_model.gap_extend),
                   std::next(std::begin(extend_penalties), hmm::min_flank_pad()),
                   [] (const auto& probs) noexcept { return probability_to_penalty(probs[1]); });
}

auto recalculate_log_probability(const CigarString& cigar, const double snv_probability,
                                 const IndelMutationModel::ContextIndelModel& indel_model)
{
    const auto snv_log_probability = std::log(snv_probability);
    double result {0};
    std::size_t pos {0};
    for (const auto& op : cigar) {
        if (op.flag() == CigarOperation::Flag::substitution) {
            result += op.size() * snv_log_probability;
            pos += op.size();
        } else if (is_indel(op)) {
            result += std::log(calculate_indel_probability(indel_model, pos, op.size()));
            if (is_deletion(op)) pos += op.size();
        } else {
            pos += op.size();
        }
    }
    return result;
}

} // namespace

DeNovoModel::LocalIndelModel DeNovoModel::generate_local_indel_model(const Haplotype& given) const
{
    LocalIndelModel result {};
    result.indel = indel_model_.evaluate(given);
    const auto num_bases = sequence_size(given);
    assert(result.indel.gap_open.size() == num_bases);
    result.open.resize(num_bases + 2 * hmm::min_flank_pad(), snv_penalty_);
    result.extend.resize(num_bases + 2 * hmm::min_flank_pad(), snv_penalty_);
    set_penalties(result.indel, result.open, result.extend);
    return result;
}

void DeNovoModel::set_local_indel_model(const unsigned given) const
{
    assert(given < gap_model_index_cache_.size());
    auto& cached_result = gap_model_index_cache_[given];
    if (!cached_result) {
        cached_result = generate_local_indel_model(haplotypes_[given]);
    }
    local_indel_model_ = std::addressof(*cached_result);
}

hmm::VariableGapExtendMutationModel DeNovoModel::make_hmm_model_from_cache() const
{
    return {snv_penalty_, local_indel_model_->open, local_indel_model_->extend};
}

void DeNovoModel::align_with_hmm(const Haplotype& target, const Haplotype& given) const
{
    pad_given(target, given, padded_given_);
    local_indel_model_->open.resize(padded_given_.size(), snv_penalty_);
    local_indel_model_->extend.resize(padded_given_.size(), snv_penalty_);
    const auto hmm_model = make_hmm_model_from_cache();
    hmm::align(target.sequence(), padded_given_, hmm_model, alignment_);
}

DeNovoModel::LogProbability
DeNovoModel::evaluate_uncached(const Haplotype& target, const Haplotype& given, const bool gap_penalties_cached) const
{
    if (!gap_penalties_cached) {
        tmp_indel_model_ = generate_local_indel_model(given);
        local_indel_model_ = std::addressof(tmp_indel_model_);
    }
    LogProbability result;
    if (can_align_with_hmm(target, given)) {
        align_with_hmm(target, given);
        result = recalculate_log_probability(alignment_.cigar, params_.snv_mutation_rate, local_indel_model_->indel);
    } else {
        result = calculate_approx_log_probability(target, given, params_.snv_mutation_rate, local_indel_model_->indel);
    }
    return min_ln_probability_ ? std::max(result, *min_ln_probability_) : result;
}

DeNovoModel::LogProbability
DeNovoModel::evaluate_uncached(const unsigned target_idx, const unsigned given_idx) const
{
    set_local_indel_model(given_idx);
    return evaluate_uncached(haplotypes_[target_idx], haplotypes_[given_idx], true);
}

DeNovoModel::LogProbability
DeNovoModel::evaluate_basic_cache(const Haplotype& target, const Haplotype& given) const
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

DeNovoModel::LogProbability
DeNovoModel::evaluate_address_cache(const Haplotype& target, const Haplotype& given) const
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
