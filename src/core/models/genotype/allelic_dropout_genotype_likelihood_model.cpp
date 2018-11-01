// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "allelic_dropout_genotype_likelihood_model.hpp"

#include <utility>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <array>
#include <limits>
#include <cassert>

#include "utils/maths.hpp"
#include "utils/concat.hpp"

namespace octopus { namespace model {

// public methods

AllelicDropoutGenotypeLikelihoodModel::Parameters
AllelicDropoutGenotypeLikelihoodModel::parameters() const
{
    return params_;
}

const HaplotypeLikelihoodArray& AllelicDropoutGenotypeLikelihoodModel::cache() const noexcept
{
    return constant_mixture_model_.cache();
}

void AllelicDropoutGenotypeLikelihoodModel::prime(const std::vector<Haplotype>& haplotypes)
{
    constant_mixture_model_.prime(haplotypes);
    variable_mixture_model_.prime(haplotypes);
}

void AllelicDropoutGenotypeLikelihoodModel::unprime() noexcept
{
    constant_mixture_model_.unprime();
    variable_mixture_model_.unprime();
}

bool AllelicDropoutGenotypeLikelihoodModel::is_primed() const noexcept
{
    assert(variable_mixture_model_.is_primed() == constant_mixture_model_.is_primed());
    return variable_mixture_model_.is_primed();
}

namespace {

template <typename T = AllelicDropoutGenotypeLikelihoodModel::LogProbability>
static constexpr auto ln(const unsigned n)
{
    constexpr std::array<T, 11> lnLookup {std::numeric_limits<T>::infinity(),
                                          0.0,
                                          0.693147180559945309417232121458176568075500134360255254120,
                                          1.098612288668109691395245236922525704647490557822749451734,
                                          1.386294361119890618834464242916353136151000268720510508241,
                                          1.609437912434100374600759333226187639525601354268517721912,
                                          1.791759469228055000812477358380702272722990692183004705855,
                                          1.945910149055313305105352743443179729637084729581861188459,
                                          2.079441541679835928251696364374529704226500403080765762362,
                                          2.197224577336219382790490473845051409294981115645498903469,
                                          2.302585092994045684017991454684364207601101488628772976033};
    return lnLookup[n];
}

template <typename G>
AllelicDropoutGenotypeLikelihoodModel::LogProbability
evaluate_dropout(const G& genotype, const unsigned ploidy, VariableMixtureGenotypeLikelihoodModel& model)
{
    VariableMixtureGenotypeLikelihoodModel::MixtureVector mixtures(ploidy, 1.0 / (ploidy - 1));
    auto result = -ln(ploidy);
    for (unsigned i {0}; i < ploidy; ++i) {
        mixtures[i] = 0;
        model.set_mixtures(mixtures);
        result += model.evaluate(genotype);
        mixtures[i] = 1.0 / (ploidy - 1);
    }
    return result;
}

} // namespace

AllelicDropoutGenotypeLikelihoodModel::LogProbability
AllelicDropoutGenotypeLikelihoodModel::evaluate(const Genotype<Haplotype>& genotype) const
{
    if (genotype.ploidy() < 2) return constant_mixture_model_.evaluate(genotype);
    return maths::log_sum_exp(std::log(dropout_probability()) + evaluate_dropout(genotype, genotype.ploidy(), variable_mixture_model_),
                              std::log(1 - dropout_probability()) + constant_mixture_model_.evaluate(genotype));
}

AllelicDropoutGenotypeLikelihoodModel::LogProbability
AllelicDropoutGenotypeLikelihoodModel::evaluate(const GenotypeIndex& genotype) const
{
    if (genotype.size() < 2) return constant_mixture_model_.evaluate(genotype);
    return maths::log_sum_exp(std::log(dropout_probability()) + evaluate_dropout(genotype, genotype.size(), variable_mixture_model_),
                              std::log(1 - dropout_probability()) + constant_mixture_model_.evaluate(genotype));
}

AllelicDropoutGenotypeLikelihoodModel::LogProbability
AllelicDropoutGenotypeLikelihoodModel::evaluate(const CancerGenotype<Haplotype>& genotype) const
{
    assert(genotype.ploidy() > 1);
    return maths::log_sum_exp(std::log(dropout_probability()) + evaluate_dropout(genotype, genotype.ploidy(), variable_mixture_model_),
                              std::log(1 - dropout_probability()) + constant_mixture_model_.evaluate(demote(genotype)));
}

AllelicDropoutGenotypeLikelihoodModel::LogProbability
AllelicDropoutGenotypeLikelihoodModel::evaluate(const CancerGenotypeIndex& genotype) const
{
    const auto merged = concat(genotype.germline, genotype.somatic);
    assert(merged.size() > 1);
    return maths::log_sum_exp(std::log(dropout_probability()) + evaluate_dropout(genotype, merged.size(), variable_mixture_model_),
                              std::log(1 - dropout_probability()) + constant_mixture_model_.evaluate(merged));
}

// private methods

AllelicDropoutGenotypeLikelihoodModel::LogProbability
AllelicDropoutGenotypeLikelihoodModel::dropout_probability() const noexcept
{
    return params_.dropout_rate;
}

} // namespace model
} // namespace octopus
