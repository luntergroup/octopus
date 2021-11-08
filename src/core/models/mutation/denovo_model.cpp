// Copyright (c) 2015-2021 Daniel Cooke
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
    return std::max(std::min(std::round(std::log(probability) / -maths::constants::ln10Div10<>), max_penalty), 1.0);
}

} // namespace

DeNovoModel::DeNovoModel(Parameters parameters, std::size_t num_haplotypes_hint, CachingStrategy caching)
: params_ {parameters}
, snv_log_prior_ {std::log(params_.snv_prior)}
, pad_penalty_ {60}
, snv_penalty_ {probability_to_penalty(params_.snv_prior)}
, indel_model_ {{params_.indel_prior}}
, min_ln_probability_ {}
, num_haplotypes_hint_ {num_haplotypes_hint}
, caching_ {caching}
, alignment_ {}
, tmp_indel_model_ {}
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

void DeNovoModel::prime(const MappableBlock<Haplotype>& haplotypes)
{
    if (is_primed()) throw std::runtime_error {"DeNovoModel: already primed"};
    constexpr std::size_t max_unguardered {50};
    gap_model_index_cache_.resize(haplotypes.size());
    if (haplotypes.size() <= max_unguardered) {
        unguarded_index_cache_.assign(haplotypes.size(), std::vector<LogProbability>(haplotypes.size(), 0));
        const auto indexed_haplotypes = index(haplotypes);
        for (std::size_t given_idx {0}; given_idx < haplotypes.size(); ++given_idx) {
            for (std::size_t target_idx {0}; target_idx < haplotypes.size(); ++target_idx) {
                if (target_idx != given_idx) {
                    unguarded_index_cache_[target_idx][given_idx] = evaluate_uncached(indexed_haplotypes[target_idx], indexed_haplotypes[given_idx]);
                }
            }
            gap_model_index_cache_[given_idx] = boost::none; // clear the cache to reclaim memory
        }
        use_unguarded_ = true;
    } else {
        guarded_index_cache_.assign(haplotypes.size(), std::vector<boost::optional<LogProbability>>(haplotypes.size()));
    }
}

void DeNovoModel::unprime() noexcept
{
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

DeNovoModel::LogProbability DeNovoModel::evaluate(const IndexedHaplotype<>& target, const IndexedHaplotype<>& given) const
{
    if (use_unguarded_) {
        return unguarded_index_cache_[index_of(target)][index_of(given)];
    } else {
        auto& result = guarded_index_cache_[index_of(target)][index_of(given)];
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
               const std::size_t flank_pad, std::string& result)
{
    const auto required_size = std::max(target.size(), given.size()) + 2 * flank_pad;
    result.resize(required_size);
    auto itr = std::fill_n(std::begin(result), flank_pad, 'N');
    itr = std::copy(std::cbegin(given), std::cend(given), itr);
    std::fill(itr, std::end(result), 'N');
}

void pad_given(const Haplotype& target, const Haplotype& given, std::size_t flank_pad, std::string& result)
{
    pad_given(target.sequence(), given.sequence(), flank_pad, result);
}

auto sequence_length_distance(const Haplotype& lhs, const Haplotype& rhs) noexcept
{
    const auto p = std::minmax({sequence_size(lhs), sequence_size(rhs)});
    return p.second - p.first;
}

DeNovoModel::LogProbability
calculate_log_probability(const Variant& variant, const Haplotype& context,
                          const DeNovoModel::LogProbability snv_log_prior, 
                          const IndelMutationModel::ContextIndelModel& indel_model)
{
    if (is_indel(variant)) {
        assert(contains(context, variant));
        const auto offset = std::min(static_cast<std::size_t>(begin_distance(context, variant)), indel_model.gap_open.size() - 1);
        if (is_simple_indel(variant)) {
            return std::log(calculate_indel_probability(indel_model, offset, indel_size(variant)));
        } else {
            return std::log(calculate_indel_probability(indel_model, offset, ref_sequence_size(variant) + alt_sequence_size(variant)));
        }
    } else {
        return snv_log_prior;
    }
}

auto calculate_approx_log_probability(const Haplotype& target, const Haplotype& given,
                                      const DeNovoModel::LogProbability snv_log_prior, 
                                      const IndelMutationModel::ContextIndelModel& indel_model)
{
    DeNovoModel::LogProbability result {0};
    const auto variants = target.difference(given);
    for (const auto& variant : variants) {
        result += calculate_log_probability(variant, given, snv_log_prior, indel_model);
    }
    return result;
}

void set_penalties(const IndelMutationModel::ContextIndelModel& indel_model,
                   const std::size_t flank_pad,
                   hmm::PenaltyVector& open_penalties,
                   hmm::PenaltyVector& extend_penalties)
{
    assert(indel_model.gap_open.size() + 2 *flank_pad == open_penalties.size());
    std::transform(std::cbegin(indel_model.gap_open), std::cend(indel_model.gap_open),
                   std::next(std::begin(open_penalties), flank_pad),
                   probability_to_penalty);
    assert(indel_model.gap_extend.size() + 2 * flank_pad == extend_penalties.size());
    std::transform(std::cbegin(indel_model.gap_extend), std::cend(indel_model.gap_extend),
                   std::next(std::begin(extend_penalties), flank_pad),
                   [] (const auto& probs) noexcept { return probability_to_penalty(probs[1]); });
}

auto recalculate_log_probability(const CigarString& alignment, 
                                 const DeNovoModel::LogProbability snv_log_prior,
                                 const IndelMutationModel::ContextIndelModel& indel_model)
{
    assert(reference_size(alignment) <= indel_model.gap_open.size());
    DeNovoModel::LogProbability result {0};
    std::size_t pos {0};
    for (const auto& op : alignment) {
        if (op.flag() == CigarOperation::Flag::substitution) {
            result += op.size() * snv_log_prior;
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
    result.open.resize(num_bases + 2 * hmm_.band_size(), pad_penalty_);
    result.extend.resize(num_bases + 2 * hmm_.band_size(), pad_penalty_);
    set_penalties(result.indel, hmm_.band_size(), result.open, result.extend);
    return result;
}

void DeNovoModel::set_local_indel_model(const IndexedHaplotype<>& given) const
{
    assert(index_of(given) < gap_model_index_cache_.size());
    auto& cached_result = gap_model_index_cache_[index_of(given)];
    if (!cached_result) {
        cached_result = generate_local_indel_model(given.haplotype());
    }
    local_indel_model_ = std::addressof(*cached_result);
}

DeNovoModel::HMM::ParameterType DeNovoModel::make_hmm_parameters() const noexcept
{
    return {local_indel_model_->open, local_indel_model_->extend, {}, {}, snv_penalty_};
}

bool DeNovoModel::can_try_align_with_hmm(const Haplotype& target, const Haplotype& given) const noexcept
{
    return sequence_length_distance(target, given) < static_cast<unsigned>(hmm_.band_size());
}

void DeNovoModel::align_with_hmm(const Haplotype& target, const Haplotype& given) const
{
    pad_given(target, given, hmm_.band_size(), padded_given_);
    local_indel_model_->open.resize(padded_given_.size(), pad_penalty_);
    local_indel_model_->extend.resize(padded_given_.size(), pad_penalty_);
    const auto hmm_params = make_hmm_parameters();
    hmm_.set(hmm_params);
    hmm_.align(target.sequence(), padded_given_, alignment_);
}

bool is_valid_alignment(const hmm::Alignment& alignment, std::size_t flank_pad) noexcept
{
    return alignment.target_offset == flank_pad;
}

DeNovoModel::LogProbability
DeNovoModel::evaluate_uncached(const Haplotype& target, const Haplotype& given, const bool gap_penalties_cached) const
{
    if (!gap_penalties_cached) {
        tmp_indel_model_ = generate_local_indel_model(given);
        local_indel_model_ = std::addressof(tmp_indel_model_);
    }
    LogProbability result;
    if (can_try_align_with_hmm(target, given)) {
        try {
            align_with_hmm(target, given);
            if (is_valid_alignment(alignment_, hmm_.band_size())) {
                result = recalculate_log_probability(alignment_.cigar, snv_log_prior_, local_indel_model_->indel);
            } else {
                result = calculate_approx_log_probability(target, given, snv_log_prior_, local_indel_model_->indel);
            }
        } catch (const hmm::HMMOverflow&) {
            result = calculate_approx_log_probability(target, given, snv_log_prior_, local_indel_model_->indel);
        }
    } else {
        result = calculate_approx_log_probability(target, given, snv_log_prior_, local_indel_model_->indel);
    }
    return min_ln_probability_ ? std::max(result, *min_ln_probability_) : result;
}

DeNovoModel::LogProbability
DeNovoModel::evaluate_uncached(const IndexedHaplotype<>& target, const IndexedHaplotype<>& given) const
{
    set_local_indel_model(given);
    return evaluate_uncached(target.haplotype(), given.haplotype(), true);
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
