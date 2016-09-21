// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "denovo_model.hpp"

#include <iterator>
#include <algorithm>

#include "basics/phred.hpp"
#include "../pairhmm/pair_hmm.hpp"

namespace octopus {

DeNovoModel::DeNovoModel(Parameters params)
: params_ {params}
{}

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
    const auto p = static_cast<std::int8_t>(probability_to_phred(params_.mutation_rate).score());
    hmm::BasicMutationModel model {p, p, p};
    return hmm::evaluate(target.sequence(), pad(given.sequence(), sequence_size(target)), model);
}

} // namespace octopus
