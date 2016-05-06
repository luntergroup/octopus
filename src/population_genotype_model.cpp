//
//  population_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "population_genotype_model.hpp"

#include <utility>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>

#include "fixed_ploidy_genotype_likelihood_model.hpp"
#include "maths.hpp"
#include "logging.hpp"

namespace Octopus
{
namespace GenotypeModel
{

Population::Population(const CoalescentModel& genotype_prior_model)
:
genotype_prior_model_ {genotype_prior_model}
{}

Population::Latents::Latents(SampleGenotypeProbabilityVector&& genotype_probabilities)
:
genotype_probabilities {std::move(genotype_probabilities)}
{}

Population::InferredLatents::InferredLatents(Latents&& posteriors, double log_evidence)
:
posteriors {std::move(posteriors)},
log_evidence {log_evidence}
{}

namespace debug
{
    template <typename S>
    void print_genotype_likelihoods(S&& stream, const std::vector<Genotype<Haplotype>>& genotypes,
                                    const std::unordered_map<SampleIdType, std::vector<double>>& likelihoods,
                                    std::size_t n = 5);
    void print_genotype_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                    const std::unordered_map<SampleIdType, std::vector<double>>& likelihoods,
                                    std::size_t n = 5);
}

namespace
{
    auto calculate_genotype_likelihoods(const std::vector<SampleIdType>& samples,
                                        const std::vector<Genotype<Haplotype>>& genotypes,
                                        const HaplotypeLikelihoodCache& haplotype_likelihoods)
    {
        std::vector<std::vector<double>> result(samples.size(), std::vector<double>(genotypes.size()));
        
        const auto ploidy = genotypes.front().ploidy();
        
        FixedPloidyGenotypeLikelihoodModel model {ploidy, haplotype_likelihoods};
        
        auto it = std::begin(result);
        for (const auto& sample : samples) {
            std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(*it++),
                           [&sample, &model] (const auto& genotype) {
                               return model.log_likelihood(sample, genotype);
                           });
        }
        
        return result;
    }
} // namespace

Population::InferredLatents
Population::infer_latents(const std::vector<SampleIdType>& samples, const GenotypeVector& genotypes,
                          const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!genotypes.empty());
    
    auto genotype_likelihoods = calculate_genotype_likelihoods(samples, genotypes, haplotype_likelihoods);
    
    std::vector<double> posteriors {};
    
    for (unsigned i = 0; i < genotypes.size(); ++i) {
        for (unsigned j = 0; j < genotypes.size(); ++j) {
            std::vector<Haplotype> haplotypes {};
            for (const auto& h : genotypes[i]) haplotypes.push_back(h);
            for (const auto& h : genotypes[j]) haplotypes.push_back(h);
            posteriors.emplace_back(genotype_prior_model_.get().evaluate(haplotypes)
                                    + genotype_likelihoods[0][i] + genotype_likelihoods[1][j]);
        }
    }
    
    Maths::normalise_exp(posteriors);
    
    Latents::SampleGenotypeProbabilityVector genotype_posteriors(samples.size(), Latents::GenotypeProbabilityVector(genotypes.size()));
    
    for (unsigned i = 0; i < genotypes.size(); ++i) {
        double p {0};
        for (unsigned j = 0; j < genotypes.size(); ++j) {
            p += posteriors[genotypes.size() * i + j];
        }
        genotype_posteriors[0][i] = p;
    }
    
    for (unsigned i = 0; i < genotypes.size(); ++i) {
        double p {0};
        for (unsigned j = 0; j < genotypes.size(); ++j) {
            p += posteriors[i + genotypes.size() * j];
        }
        genotype_posteriors[1][i] = p;
    }
    
    InferredLatents result {};
    
    return result;
}

namespace debug
{
    template <typename S>
    void print_genotype_likelihoods(S&& stream, const std::vector<Genotype<Haplotype>>& genotypes,
                                    const std::unordered_map<SampleIdType, std::vector<double>>& likelihoods,
                                    std::size_t n)
    {
//        assert(genotypes.size() == likelihoods.size());
//        
//        const auto m = std::min(n, genotypes.size());
//        
//        if (m == genotypes.size()) {
//            stream << "Printing all genotype likelihoods " << '\n';
//        } else {
//            stream << "Printing top " << m << " genotype likelihoods " << '\n';
//        }
//        
//        using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;
//        
//        std::vector<std::pair<GenotypeReference, double>> v {};
//        v.reserve(genotypes.size());
//        
//        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(likelihoods),
//                       std::back_inserter(v), [] (const auto& g, const auto& p) {
//                           return std::make_pair(std::cref(g), p);
//                       });
//        
//        const auto mth = std::next(std::begin(v), m);
//        
//        std::partial_sort(std::begin(v), mth, std::end(v),
//                          [] (const auto& lhs, const auto& rhs) {
//                              return lhs.second > rhs.second;
//                          });
//        
//        std::for_each(std::begin(v), mth,
//                      [&] (const auto& p) {
//                          ::debug::print_variant_alleles(stream, p.first);
//                          stream << " " << p.second << '\n';
//                      });
    }
    
    void print_genotype_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                    const std::unordered_map<SampleIdType, std::vector<double>>& likelihoods,
                                    std::size_t n)
    {
        print_genotype_likelihoods(std::cout, genotypes, likelihoods, n);
    }
} // namespace debug
} // namespace GenotypeModel
} // namespace Octopus

