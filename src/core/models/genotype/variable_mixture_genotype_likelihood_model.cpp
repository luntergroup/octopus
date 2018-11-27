// Copyright (c) 2015-2018 Daniel Cooke
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
, indexed_likelihoods_ {}
, likelihood_refs_ {}
, buffer_ {}
{}

VariableMixtureGenotypeLikelihoodModel::VariableMixtureGenotypeLikelihoodModel(const HaplotypeLikelihoodArray& likelihoods,
                                                                         const std::vector<Haplotype>& haplotypes)
: VariableMixtureGenotypeLikelihoodModel {likelihoods}
{
    this->prime(haplotypes);
}

VariableMixtureGenotypeLikelihoodModel::VariableMixtureGenotypeLikelihoodModel(const HaplotypeLikelihoodArray& likelihoods, MixtureVector mixtures)
: VariableMixtureGenotypeLikelihoodModel {likelihoods}
{
    this->set_mixtures(std::move(mixtures));
}

VariableMixtureGenotypeLikelihoodModel::VariableMixtureGenotypeLikelihoodModel(const HaplotypeLikelihoodArray& likelihoods,
                                                                         MixtureVector mixtures,
                                                                         const std::vector<Haplotype>& haplotypes)
: VariableMixtureGenotypeLikelihoodModel {likelihoods, std::move(mixtures)}
{
    this->prime(haplotypes);
}

const HaplotypeLikelihoodArray& VariableMixtureGenotypeLikelihoodModel::cache() const noexcept
{
    return likelihoods_;
}

const VariableMixtureGenotypeLikelihoodModel::MixtureVector& VariableMixtureGenotypeLikelihoodModel::mixtures() const noexcept
{
    return mixtures_;
}

void VariableMixtureGenotypeLikelihoodModel::prime(const std::vector<Haplotype>& haplotypes)
{
    assert(likelihoods_.is_primed());
    indexed_likelihoods_.reserve(haplotypes.size());
    std::transform(std::cbegin(haplotypes), std::cend(haplotypes), std::back_inserter(indexed_likelihoods_),
                   [this] (const auto& haplotype) -> const HaplotypeLikelihoodArray::LikelihoodVector& {
                       return likelihoods_[haplotype]; });
}

void VariableMixtureGenotypeLikelihoodModel::unprime() noexcept
{
    indexed_likelihoods_.clear();
    indexed_likelihoods_.shrink_to_fit();
}

bool VariableMixtureGenotypeLikelihoodModel::is_primed() const noexcept
{
    return !indexed_likelihoods_.empty();
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
VariableMixtureGenotypeLikelihoodModel::evaluate(const GenotypeIndex& genotype) const
{
    assert(is_primed());
    assert(genotype.size() == mixtures_.size());
    assert(buffer_.size() == mixtures_.size());
    LogProbability result {0};
    const auto num_reads = indexed_likelihoods_.front().get().size();
    for (std::size_t read_idx {0}; read_idx < num_reads; ++read_idx) {
        std::transform(std::cbegin(genotype), std::cend(genotype), std::cbegin(log_mixtures_), std::begin(buffer_),
                       [this, read_idx] (auto haplotype_idx, auto log_mixture) noexcept {
                           return log_mixture + indexed_likelihoods_[haplotype_idx].get()[read_idx]; });
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
VariableMixtureGenotypeLikelihoodModel::evaluate(const CancerGenotypeIndex& genotype) const
{
    assert(is_primed());
    assert((genotype.germline.size() + genotype.somatic.size()) == mixtures_.size());
    assert(buffer_.size() == mixtures_.size());
    LogProbability result {0};
    const auto num_reads = indexed_likelihoods_.front().get().size();
    for (std::size_t read_idx {0}; read_idx < num_reads; ++read_idx) {
        auto buffer_itr = std::transform(std::cbegin(genotype.germline), std::cend(genotype.germline),
                                         std::cbegin(log_mixtures_), std::begin(buffer_),
                                         [=] (auto haplotype_idx, auto log_mixture) noexcept {
                                             return log_mixture + indexed_likelihoods_[haplotype_idx].get()[read_idx]; });
        std::transform(std::cbegin(genotype.somatic), std::cend(genotype.somatic),
                       std::next(std::cbegin(log_mixtures_), genotype.germline.size()), buffer_itr,
                       [this, read_idx] (auto haplotype_idx, auto log_mixture) noexcept {
                           return log_mixture + indexed_likelihoods_[haplotype_idx].get()[read_idx]; });
        result += maths::log_sum_exp(buffer_);
    }
    return result;
}

} // namespace model
} // namespace octopus
