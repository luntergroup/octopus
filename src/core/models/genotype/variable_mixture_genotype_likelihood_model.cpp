// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "variable_mixture_genotype_likelihood_model.hpp"

#include <utility>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <array>
#include <limits>
#include <cassert>

#include "utils/maths.hpp"

namespace octopus { namespace model {

VariableMixtureGenotypeLikelihoodModel::VariableMixtureGenotypeLikelihoodModel(const HaplotypeLikelihoodArray& likelihoods)
: likelihoods_ {likelihoods}
, mixtures_ {}
, likelihood_refs_ {}
, buffer_ {}
{}

VariableMixtureGenotypeLikelihoodModel::VariableMixtureGenotypeLikelihoodModel(const HaplotypeLikelihoodArray& likelihoods, MixtureVector mixtures)
: VariableMixtureGenotypeLikelihoodModel {likelihoods}
{
    this->set_mixtures(std::move(mixtures));
}

const HaplotypeLikelihoodArray& VariableMixtureGenotypeLikelihoodModel::cache() const noexcept
{
    return likelihoods_;
}

const VariableMixtureGenotypeLikelihoodModel::MixtureVector& VariableMixtureGenotypeLikelihoodModel::mixtures() const noexcept
{
    return mixtures_;
}

void VariableMixtureGenotypeLikelihoodModel::set_mixtures(MixtureVector mixtures)
{
    mixtures_ = std::move(mixtures);
    log_mixtures_ = mixtures_;
    maths::log_each(log_mixtures_);
    likelihood_refs_.reserve(mixtures_.size());
    buffer_.resize(mixtures_.size());
}

VariableMixtureGenotypeLikelihoodModel::LogProbability
VariableMixtureGenotypeLikelihoodModel::evaluate(const Genotype<Haplotype>& genotype) const
{
    assert(genotype.ploidy() == mixtures_.size());
    assert(buffer_.size() == mixtures_.size());
    likelihood_refs_.clear();
    std::transform(std::cbegin(genotype), std::cend(genotype), std::back_inserter(likelihood_refs_),
                   [this] (const auto& haplotype) -> const HaplotypeLikelihoodArray::LikelihoodVector& {
                       return likelihoods_[haplotype]; });
    LogProbability result {0};
    const auto num_reads = likelihood_refs_.front().get().size();
    for (std::size_t read_idx {0}; read_idx < num_reads; ++read_idx) {
        std::transform(std::cbegin(likelihood_refs_), std::cend(likelihood_refs_),
                       std::cbegin(log_mixtures_), std::begin(buffer_),
                       [read_idx] (const auto& likelihoods, auto log_mixture) noexcept {
                            return log_mixture + likelihoods.get()[read_idx]; });
        result += maths::log_sum_exp(buffer_);
    }
    return result;
}

VariableMixtureGenotypeLikelihoodModel::LogProbability
VariableMixtureGenotypeLikelihoodModel::evaluate(const Genotype<IndexedHaplotype<>>& genotype) const
{
    assert(genotype.ploidy() == mixtures_.size());
    LogProbability result {0};
    const auto num_likelihoods = likelihoods_.num_likelihoods();
    for (std::size_t read_idx {0}; read_idx < num_likelihoods; ++read_idx) {
        std::transform(std::cbegin(genotype), std::cend(genotype), std::cbegin(log_mixtures_), std::begin(buffer_),
                       [this, read_idx] (const auto& haplotype, auto log_mixture) noexcept {
                           return log_mixture + likelihoods_[haplotype][read_idx]; });
        result += maths::log_sum_exp(buffer_);
    }
    return result;
}

VariableMixtureGenotypeLikelihoodModel::LogProbability
VariableMixtureGenotypeLikelihoodModel::evaluate(const CancerGenotype<Haplotype>& genotype) const
{
    assert(genotype.ploidy() == mixtures_.size());
    assert(buffer_.size() == mixtures_.size());
    likelihood_refs_.clear();
    std::transform(std::cbegin(genotype.germline()), std::cend(genotype.germline()), std::back_inserter(likelihood_refs_),
                  [this] (const auto& haplotype) -> const HaplotypeLikelihoodArray::LikelihoodVector& { return likelihoods_[haplotype]; });
    std::transform(std::cbegin(genotype.somatic()), std::cend(genotype.somatic()), std::back_inserter(likelihood_refs_),
                   [this] (const auto& haplotype) -> const HaplotypeLikelihoodArray::LikelihoodVector& { return likelihoods_[haplotype]; });
    LogProbability result {0};
    const auto num_reads = likelihood_refs_.front().get().size();
    for (std::size_t read_idx {0}; read_idx < num_reads; ++read_idx) {
        std::transform(std::cbegin(likelihood_refs_), std::cend(likelihood_refs_),
                       std::cbegin(log_mixtures_), std::begin(buffer_),
                       [read_idx] (const auto& likelihoods, auto log_mixture) noexcept {
                           return log_mixture + likelihoods.get()[read_idx]; });
        result += maths::log_sum_exp(buffer_);
    }
    return result;
}

VariableMixtureGenotypeLikelihoodModel::LogProbability
VariableMixtureGenotypeLikelihoodModel::evaluate(const CancerGenotype<IndexedHaplotype<>>& genotype) const
{
    assert((genotype.germline_ploidy() + genotype.somatic_ploidy()) == mixtures_.size());
    assert(buffer_.size() == mixtures_.size());
    LogProbability result {0};
    const auto num_likelihoods = likelihoods_.num_likelihoods();
    for (std::size_t read_idx {0}; read_idx < num_likelihoods; ++read_idx) {
        auto buffer_itr = std::transform(std::cbegin(genotype.germline()), std::cend(genotype.germline()),
                                         std::cbegin(log_mixtures_), std::begin(buffer_),
                                         [=] (const auto& haplotype, auto log_mixture) noexcept {
                                             return log_mixture + likelihoods_[haplotype][read_idx]; });
        std::transform(std::cbegin(genotype.somatic()), std::cend(genotype.somatic()),
                       std::next(std::cbegin(log_mixtures_), genotype.germline_ploidy()), buffer_itr,
                       [this, read_idx] (const auto& haplotype, auto log_mixture) noexcept {
                           return log_mixture + likelihoods_[haplotype][read_idx]; });
        result += maths::log_sum_exp(buffer_);
    }
    return result;
}

} // namespace model
} // namespace octopus
