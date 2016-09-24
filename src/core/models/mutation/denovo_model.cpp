// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "denovo_model.hpp"

#include <iterator>
#include <algorithm>
#include <cstdint>
#include <cassert>

#include "basics/phred.hpp"
#include "../pairhmm/pair_hmm.hpp"

namespace octopus {

DeNovoModel::DeNovoModel(Parameters parameters, std::size_t num_haplotypes_hint)
: parameters_ {parameters}
, num_haplotypes_hint_ {num_haplotypes_hint}
, cache_ {num_haplotypes_hint}
{}

auto make_hmm_model(const double denovo_mutation_rate) noexcept
{
    const auto p = static_cast<std::int8_t>(probability_to_phred(denovo_mutation_rate).score());
    return hmm::BasicMutationModel {p, p, p};
}

auto pad(const Haplotype::NucleotideSequence& given, const std::size_t target_size)
{
    auto required_pad = 2 * hmm::min_flank_pad();
    const auto given_size = given.size();
    if (target_size > given_size) {
        required_pad += target_size - given_size;
    } else if (given_size > target_size && (given_size - target_size) > required_pad) {
        return given;
    }
    Haplotype::NucleotideSequence result(given.size() + required_pad, 'N');
    std::copy(std::cbegin(given), std::cend(given),
              std::next(std::begin(result), hmm::min_flank_pad()));
    return result;
}

double DeNovoModel::evaluate(const Haplotype& target, const Haplotype& given) const
{
    const auto target_iter = cache_.find(target);
    if (target_iter != std::cend(cache_)) {
        const auto given_iter = target_iter->second.find(given);
        if (given_iter != std::cend(target_iter->second)) {
            return given_iter->second;
        }
    }
    const auto model = make_hmm_model(parameters_.mutation_rate);
    const auto result = hmm::evaluate(target.sequence(), pad(given.sequence(), sequence_size(target)), model);
    if (target_iter != std::cend(cache_)) {
        target_iter->second.emplace(given, result);
    } else {
        auto p = cache_.emplace(std::piecewise_construct,
                                std::forward_as_tuple(target),
                                std::forward_as_tuple(num_haplotypes_hint_));
        assert(p.second);
        p.first->second.emplace(given, result);
    }
    return result;
}

} // namespace octopus
