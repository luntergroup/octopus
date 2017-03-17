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

#include "basics/phred.hpp"
#include "core/types/variant.hpp"
//#include "../pairhmm/pair_hmm.hpp"

namespace octopus {

DeNovoModel::DeNovoModel(Parameters parameters, std::size_t num_haplotypes_hint, CachingStrategy caching)
: parameters_ {parameters}
, num_haplotypes_hint_ {num_haplotypes_hint}
, caching_ {caching}
, value_cache_ {}
, address_cache_ {}
{
    if (caching_ == CachingStrategy::address) {
        address_cache_.reserve(num_haplotypes_hint_ * num_haplotypes_hint_);
    } else if (caching == CachingStrategy::value) {
        value_cache_.reserve(num_haplotypes_hint_);
    }
}

//auto make_hmm_model(const double denovo_mutation_rate) noexcept
//{
//    const auto p = static_cast<std::int8_t>(probability_to_phred(denovo_mutation_rate).score());
//    return hmm::BasicMutationModel {p, p, p};
//}
//
//auto pad(const Haplotype::NucleotideSequence& given, const std::size_t target_size)
//{
//    auto required_pad = 2 * hmm::min_flank_pad();
//    const auto given_size = given.size();
//    if (target_size > given_size) {
//        required_pad += target_size - given_size;
//    } else if (given_size > target_size) {
//        const auto excess = given_size - target_size;
//        if (excess >= required_pad) {
//            return given;
//        } else {
//            required_pad -= excess;
//        }
//    }
//    Haplotype::NucleotideSequence result(given.size() + required_pad, 'N');
//    std::copy(std::cbegin(given), std::cend(given),
//              std::next(std::begin(result), hmm::min_flank_pad()));
//    return result;
//}

double DeNovoModel::evaluate(const Haplotype& target, const Haplotype& given) const
{
    switch (caching_) {
        case CachingStrategy::address: return evaluate_address_cache(target, given);
        case CachingStrategy::value: return evaluate_basic_cache(target, given);
        case CachingStrategy::none: return evaluate_uncached(target, given);
        default: return evaluate_uncached(target, given); // to prevent compiler warning
    }
}

// private methods

double DeNovoModel::evaluate_uncached(const Haplotype& target, const Haplotype& given) const
{
    // TODO: make indel errors context based
    const auto variants = target.difference(given);
    return std::accumulate(std::cbegin(variants), std::cend(variants), 0.0,
                           [this] (const double curr, const Variant& v) {
                               double penalty {std::log(parameters_.mutation_rate)};
//                               if (is_insertion(v)) {
//                                   penalty *= alt_sequence_size(v);
//                               } else if (is_deletion(v)) {
//                                   penalty *= region_size(v);
//                               }
                               return curr + penalty;
                           });
//    const auto model = make_hmm_model(parameters_.mutation_rate);
//    const auto result = hmm::evaluate(target.sequence(), pad(given.sequence(), sequence_size(target)), model);
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
    if (first_haplotype_) {
        return evaluate_flat_address_cache(target, given);
    } else {
        return evaluate_sparse_address_cache(target, given);
    }
}

double DeNovoModel::evaluate_flat_address_cache(const Haplotype& target, const Haplotype& given) const noexcept
{
    assert(first_haplotype_);
    assert(*first_haplotype_ <= std::addressof(target) && std::addressof(target) <= *last_haplotype_);
    assert(*first_haplotype_ <= std::addressof(given) && std::addressof(given) <= *last_haplotype_);
    const auto target_index = static_cast<std::size_t>(std::distance(*first_haplotype_, std::addressof(target)));
    const auto given_index  = static_cast<std::size_t>(std::distance(*first_haplotype_, std::addressof(given)));
    assert(target_index < flat_address_cache_.size());
    assert(given_index < flat_address_cache_.size());
    auto& result = flat_address_cache_[target_index][given_index];
    if (!result) {
        result = evaluate_uncached(target, given);
    }
    return *result;
}

double DeNovoModel::evaluate_sparse_address_cache(const Haplotype& target, const Haplotype& given) const
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

void DeNovoModel::do_prime(const std::vector<const Haplotype*>& addresses) const
{
    if (!addresses.empty()) {
        const auto p = std::minmax_element(std::cbegin(addresses), std::cend(addresses));
        const auto num_spaces = static_cast<std::size_t>(std::distance(*p.first, *p.second)) + 1;
        if (num_spaces <= max_flat_spaces_) {
            first_haplotype_ = *p.first;
            last_haplotype_  = *p.second;
            flat_address_cache_.assign(num_spaces, std::vector<boost::optional<double>>(num_spaces, boost::none));
        }
    }
}

void DeNovoModel::unprime() const noexcept
{
    first_haplotype_ = boost::none;
    last_haplotype_  = boost::none;
    flat_address_cache_.clear();
    flat_address_cache_.shrink_to_fit();
}

} // namespace octopus
