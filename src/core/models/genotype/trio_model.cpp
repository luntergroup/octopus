// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "trio_model.hpp"

#include <functional>
#include <cstddef>
#include <iterator>
#include <algorithm>
#include <cassert>

#include "utils/maths.hpp"
#include "germline_likelihood_model.hpp"

namespace octopus { namespace model {

TrioModel::TrioModel(const Trio& trio,
              const CoalescentModel& genotype_prior_model,
              const DeNovoModel& mutation_model,
              boost::optional<logging::DebugLogger> debug_log)
: trio_ {trio}
, genotype_prior_model_ {genotype_prior_model}
, mutation_model_ {mutation_model}
, debug_log_ {debug_log}
{}

auto compute_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                         const GermlineLikelihoodModel& model)
{
    std::vector<double> result(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                   [&model] (const auto& genotype) { return model.evaluate(genotype); });
    return result;
}

struct GenotypeRefProbabilityPair
{
    std::reference_wrapper<const Genotype<Haplotype>> genotype;
    double probability;
};

bool operator<(const GenotypeRefProbabilityPair& lhs, const GenotypeRefProbabilityPair& rhs) noexcept
{
    return lhs.probability < rhs.probability;
}
bool operator>(const GenotypeRefProbabilityPair& lhs, const GenotypeRefProbabilityPair& rhs) noexcept
{
    return lhs.probability > rhs.probability;
}

auto zip(const std::vector<Genotype<Haplotype>>& genotypes,
         const std::vector<double>& probabilities)
{
    assert(genotypes.size() == probabilities.size());
    std::vector<GenotypeRefProbabilityPair> result {};
    result.reserve(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes),
                   std::cbegin(probabilities), std::back_inserter(result),
                   [] (const auto& genotype, const auto p) {
                       return GenotypeRefProbabilityPair {genotype, p};
                   });
    return result;
}

void erase_least_likely(std::vector<GenotypeRefProbabilityPair>& zipped, const std::size_t max_keep)
{
    const auto first_erase = std::next(std::begin(zipped), std::min(zipped.size(), max_keep));
    std::partial_sort(std::begin(zipped), first_erase, std::end(zipped), std::greater<> {});
    zipped.erase(first_erase, std::end(zipped));
}

TrioModel::InferredLatents
TrioModel::evaluate(const GenotypeVector& maternal_genotypes,
                    const GenotypeVector& paternal_genotypes,
                    const GenotypeVector& child_genotypes,
                    const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!maternal_genotypes.empty());
    assert(!paternal_genotypes.empty());
    assert(!child_genotypes.empty());
    
    const GermlineLikelihoodModel likelihood_model {haplotype_likelihoods};
    const auto maternal_likelihoods = compute_likelihoods(maternal_genotypes, likelihood_model);
    const auto paternal_likelihoods = compute_likelihoods(paternal_genotypes, likelihood_model);
    
    auto maternal_sorted_likelihoods = zip(maternal_genotypes, maternal_likelihoods);
    auto paternal_sorted_likelihoods = zip(paternal_genotypes, paternal_likelihoods);
    constexpr std::size_t max_genotypes {1000};
    erase_least_likely(maternal_sorted_likelihoods, max_genotypes);
    erase_least_likely(paternal_sorted_likelihoods, max_genotypes);
    
    
    
    const auto child_likelihoods = compute_likelihoods(child_genotypes, likelihood_model);
    
    InferredLatents result {};
    
    return result;
}

} // namespace model
} // namespace octopus
