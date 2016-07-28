//
//  individual_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "individual_genotype_model.hpp"

#include <utility>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>

#include "fixed_ploidy_genotype_likelihood_model.hpp"
#include "maths.hpp"

#include "timers.hpp"

namespace octopus
{
namespace model
{
Individual::Individual(const CoalescentModel& genotype_prior_model,
                       boost::optional<Logging::DebugLogger> debug_log)
:
genotype_prior_model_ {genotype_prior_model},
debug_log_ {debug_log}
{}

Individual::Latents::Latents(GenotypeProbabilityVector&& genotype_probabilities)
:
genotype_probabilities {std::move(genotype_probabilities)}
{}

Individual::InferredLatents::InferredLatents(Latents&& posteriors, double log_evidence)
:
posteriors {std::move(posteriors)},
log_evidence {log_evidence}
{}

namespace debug
{
    template <typename S>
    void print_genotype_priors(S&& stream, const std::vector<Genotype<Haplotype>>& genotypes,
                               const std::vector<double>& priors, std::size_t n = 5);
    void print_genotype_priors(const std::vector<Genotype<Haplotype>>& genotypes,
                               const std::vector<double>& priors, std::size_t n = 5);
    template <typename S>
    void print_genotype_likelihoods(S&& stream, const std::vector<Genotype<Haplotype>>& genotypes,
                                    const std::vector<double>& likelihoods, std::size_t n = 5);
    void print_genotype_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                    const std::vector<double>& likelihoods, std::size_t n = 5);
}

Individual::InferredLatents
Individual::infer_latents(const SampleName& sample,
                          const std::vector<Genotype<Haplotype>>& genotypes,
                          const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    using std::cbegin; using std::cend; using std::begin;
    
    assert(!genotypes.empty());
    
    const auto ploidy = genotypes.front().ploidy();
    
    const FixedPloidyGenotypeLikelihoodModel likelihood_model {ploidy, haplotype_likelihoods};
    
    std::vector<double> result(genotypes.size());
    
    std::transform(cbegin(genotypes), cend(genotypes), begin(result),
                   [&sample, &likelihood_model] (const auto& genotype) {
                       return likelihood_model.log_likelihood(sample, genotype);
                   });
    
    if (debug_log_) debug::print_genotype_likelihoods(stream(*debug_log_), genotypes, result);
    
    std::transform(cbegin(genotypes), cend(genotypes), cbegin(result), begin(result),
                   [this] (const auto& genotype, const auto likelihood) {
                       return genotype_prior_model_.evaluate(genotype) + likelihood;
                   });
    
    auto log_evidence = Maths::normalise_exp(result);
    
    return InferredLatents {Latents {std::move(result)}, log_evidence};
}

namespace debug
{
    template <typename S>
    void print_genotype_priors(S&& stream, const std::vector<Genotype<Haplotype>>& genotypes,
                               const std::vector<double>& priors, const std::size_t n)
    {
        assert(genotypes.size() == priors.size());
        
        const auto m = std::min(n, genotypes.size());
        
        if (m == genotypes.size()) {
            stream << "Printing all genotype priors " << '\n';
        } else {
            stream << "Printing top " << m << " genotype priors " << '\n';
        }
        
        using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;
        
        std::vector<std::pair<GenotypeReference, double>> v {};
        v.reserve(genotypes.size());
        
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(priors),
                       std::back_inserter(v), [] (const auto& g, const auto& p) {
                           return std::make_pair(std::cref(g), p);
                       });
        
        const auto mth = std::next(std::begin(v), m);
        
        std::partial_sort(std::begin(v), mth, std::end(v),
                          [] (const auto& lhs, const auto& rhs) {
                              return lhs.second > rhs.second;
                          });
        
        std::for_each(std::begin(v), mth,
                      [&] (const auto& p) {
                          ::debug::print_variant_alleles(stream, p.first);
                          stream << " " << p.second << '\n';
                      });
    }
    
    void print_genotype_priors(const std::vector<Genotype<Haplotype>>& genotypes,
                               const std::vector<double>& priors, const std::size_t n)
    {
        print_genotype_priors(std::cout, genotypes, priors, n);
    }
    
    template <typename S>
    void print_genotype_likelihoods(S&& stream, const std::vector<Genotype<Haplotype>>& genotypes,
                                    const std::vector<double>& likelihoods, std::size_t n)
    {
        assert(genotypes.size() == likelihoods.size());
        
        const auto m = std::min(n, genotypes.size());
        
        if (m == genotypes.size()) {
            stream << "Printing all genotype likelihoods " << '\n';
        } else {
            stream << "Printing top " << m << " genotype likelihoods " << '\n';
        }
        
        using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;
        
        std::vector<std::pair<GenotypeReference, double>> v {};
        v.reserve(genotypes.size());
        
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(likelihoods),
                       std::back_inserter(v), [] (const auto& g, const auto& p) {
                           return std::make_pair(std::cref(g), p);
                       });
        
        const auto mth = std::next(std::begin(v), m);
        
        std::partial_sort(std::begin(v), mth, std::end(v),
                          [] (const auto& lhs, const auto& rhs) {
                              return lhs.second > rhs.second;
                          });
        
        std::for_each(std::begin(v), mth,
                      [&] (const auto& p) {
                          ::debug::print_variant_alleles(stream, p.first);
                          stream << " " << p.second << '\n';
                      });
    }
    
    void print_genotype_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                    const std::vector<double>& likelihoods, std::size_t n)
    {
        print_genotype_likelihoods(std::cout, genotypes, likelihoods, n);
    }
} // namespace debug
} // namesapce GenotypeModel
} // namespace octopus
