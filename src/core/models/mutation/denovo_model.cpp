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

auto pad(const Haplotype::NucleotideSequence& sequence)
{
    static const auto required_pad = 2 * hmm::min_flank_pad();
    Haplotype::NucleotideSequence result(sequence.size() + required_pad, 'N');
    std::copy(std::cbegin(sequence), std::cend(sequence),
              std::next(std::begin(result), required_pad / 2));
    return result;
}

double DeNovoModel::evaluate(const Haplotype& target, const Haplotype& given) const
{
    const auto p = static_cast<std::int8_t>(probability_to_phred(params_.mutation_rate).score());
    hmm::BasicMutationModel model {p, p, p};
    return hmm::evaluate(target.sequence(), pad(given.sequence()), model);
}

} // namespace octopus
