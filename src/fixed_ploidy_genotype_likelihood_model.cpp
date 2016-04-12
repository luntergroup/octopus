//
//  fixed_ploidy_genotype_likelihood_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "fixed_ploidy_genotype_likelihood_model.hpp"

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

#include "maths.hpp"

#include <iostream> // TEST

namespace Octopus
{
namespace GenotypeModel
{
    FixedPloidyGenotypeLikelihoodModel::FixedPloidyGenotypeLikelihoodModel(unsigned ploidy,
                                                                           const HaplotypeLikelihoodCache& haplotype_likelihoods)
    :
    haplotype_likelihoods_ {haplotype_likelihoods},
    ploidy_ {ploidy},
    ln_ploidy_ {std::log(ploidy)}
    {}
    
    // ln p(read | genotype)  = ln sum {haplotype in genotype} p(read | haplotype) - ln ploidy
    // ln p(reads | genotype) = sum {read in reads} ln p(read | genotype)
    double FixedPloidyGenotypeLikelihoodModel::log_likelihood(const SampleIdType& sample,
                                                              const Genotype<Haplotype>& genotype) const
    {
        // These cases are just for optimisation
        switch (ploidy_) {
            case 1:
                return log_likelihood_haploid(sample, genotype);
            case 2:
                return log_likelihood_diploid(sample, genotype);
            case 3:
                return log_likelihood_triploid(sample, genotype);
            case 4:
                return log_likelihood_polyploid(sample, genotype);
                //return log_likelihood_tetraploid(sample, genotype);
            default:
                return log_likelihood_polyploid(sample, genotype);
        }
    }
    
    double FixedPloidyGenotypeLikelihoodModel::log_likelihood_haploid(const SampleIdType& sample,
                                                                      const Genotype<Haplotype>& genotype) const
    {
        const auto& log_likelihoods = haplotype_likelihoods_.get().log_likelihoods(sample, genotype[0]);
        return std::accumulate(std::cbegin(log_likelihoods), std::cend(log_likelihoods), 0.0);
    }
    
    double FixedPloidyGenotypeLikelihoodModel::log_likelihood_diploid(const SampleIdType& sample,
                                                                      const Genotype<Haplotype>& genotype) const
    {
        const auto& log_likelihoods1 = haplotype_likelihoods_.get().log_likelihoods(sample, genotype[0]);
        
        if (genotype.is_homozygous()) {
            return std::accumulate(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1), 0.0);
        }
        
        const auto& log_likelihoods2 = haplotype_likelihoods_.get().log_likelihoods(sample, genotype[1]);
        
        return std::inner_product(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1),
                                  std::cbegin(log_likelihoods2), 0.0, std::plus<void> {},
                                  [this] (const auto a, const auto b) -> double {
                                      return Maths::log_sum_exp(a, b) - ln_ploidy_;
                                  });
    }
    
    double FixedPloidyGenotypeLikelihoodModel::log_likelihood_triploid(const SampleIdType& sample,
                                                                       const Genotype<Haplotype>& genotype) const
    {
        const auto& log_likelihoods1 = haplotype_likelihoods_.get().log_likelihoods(sample, genotype[0]);
        
        if (genotype.is_homozygous()) {
            return std::accumulate(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1), 0.0);
        }
        
        if (genotype.zygosity() == 3) {
            const auto& log_likelihoods2 = haplotype_likelihoods_.get().log_likelihoods(sample, genotype[1]);
            const auto& log_likelihoods3 = haplotype_likelihoods_.get().log_likelihoods(sample, genotype[2]);
            return Maths::inner_product(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1),
                                        std::cbegin(log_likelihoods2), std::cbegin(log_likelihoods3),
                                        0.0, std::plus<void> {},
                                        [this] (const auto a, const auto b, const auto c) -> double {
                                            return Maths::log_sum_exp(a, b, c) - ln_ploidy_;
                                        });
        }
        
        static const double ln2 {std::log(2)};
        
        if (genotype[0] != genotype[1]) {
            const auto& log_likelihoods2 = haplotype_likelihoods_.get().log_likelihoods(sample, genotype[1]);
            return std::inner_product(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1),
                                      std::cbegin(log_likelihoods2), 0.0, std::plus<void> {},
                                      [this] (const auto a, const auto b) -> double {
                                          return Maths::log_sum_exp(a, ln2 + b) - ln_ploidy_;
                                      });
        }
        
        const auto& log_likelihoods3 = haplotype_likelihoods_.get().log_likelihoods(sample, genotype[2]);
        
        return std::inner_product(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1),
                                  std::cbegin(log_likelihoods3), 0.0, std::plus<void> {},
                                  [this] (const auto a, const auto b) -> double {
                                      return Maths::log_sum_exp(ln2 + a, b) - ln_ploidy_;
                                  });
    }
    
    double FixedPloidyGenotypeLikelihoodModel::log_likelihood_tetraploid(const SampleIdType& sample,
                                                                         const Genotype<Haplotype>& genotype) const
    {
        const auto z = genotype.zygosity();
        
        const auto& log_likelihoods1 = haplotype_likelihoods_.get().log_likelihoods(sample, genotype[0]);
        
        if (z == 1) {
            return std::accumulate(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1), 0.0);
        }
        
        if (z == 4) {
            const auto& log_likelihoods2 = haplotype_likelihoods_.get().log_likelihoods(sample, genotype[1]);
            const auto& log_likelihoods3 = haplotype_likelihoods_.get().log_likelihoods(sample, genotype[2]);
            const auto& log_likelihoods4 = haplotype_likelihoods_.get().log_likelihoods(sample, genotype[3]);
            return Maths::inner_product(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1),
                                        std::cbegin(log_likelihoods2), std::cbegin(log_likelihoods3),
                                        std::cbegin(log_likelihoods4), 0.0, std::plus<void> {},
                                        [this] (const auto a, const auto b, const auto c, const auto d) -> double {
                                            return Maths::log_sum_exp({a, b, c, d}) - ln_ploidy_;
                                        });
        }
        
        const auto count1 = genotype.count(genotype[0]);
        
        // TODO
        
        return 0;
    }
    
    double FixedPloidyGenotypeLikelihoodModel::log_likelihood_polyploid(const SampleIdType& sample,
                                                                        const Genotype<Haplotype>& genotype) const
    {
        const auto z = genotype.zygosity();
        
        const auto& log_likelihoods1 = haplotype_likelihoods_.get().log_likelihoods(sample, genotype[0]);
        
        if (z == 1) {
            return std::accumulate(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1), 0.0);
        }
        
        if (z == 2) {
            static const double lnpm1 {std::log(ploidy_ - 1)};
            
            const auto unique_haplotypes = genotype.copy_unique_ref();
            
            const auto& log_likelihoods2 = haplotype_likelihoods_.get().log_likelihoods(sample, unique_haplotypes.back());
            
            if (genotype.count(unique_haplotypes.front()) == 1) {
                return std::inner_product(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1),
                                          std::cbegin(log_likelihoods2), 0.0, std::plus<void> {},
                                          [this] (const auto a, const auto b) -> double {
                                              return Maths::log_sum_exp(a, lnpm1 + b) - ln_ploidy_;
                                          });
            }
            
            return std::inner_product(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1),
                                      std::cbegin(log_likelihoods2), 0.0, std::plus<void> {},
                                      [this] (const auto a, const auto b) -> double {
                                          return Maths::log_sum_exp(lnpm1 + a, b) - ln_ploidy_;
                                      });
        }
        
        std::vector<std::reference_wrapper<const HaplotypeLikelihoodCache::Likelihoods>> log_likelihoods {};
        log_likelihoods.reserve(ploidy_);
        
        std::transform(std::cbegin(genotype), std::cend(genotype),
                       std::back_inserter(log_likelihoods),
                       [this, &sample] (const auto& haplotype)
                            -> const HaplotypeLikelihoodCache::Likelihoods& {
                           return haplotype_likelihoods_.get().log_likelihoods(sample, haplotype);
                       });
        
        std::vector<double> tmp(ploidy_);
        
        double result {0};
        
        const auto num_likelihoods = log_likelihoods.front().get().size();
        
        for (std::size_t i {0}; i < num_likelihoods; ++i) {
            std::transform(std::cbegin(log_likelihoods), std::cend(log_likelihoods), std::begin(tmp),
                           [i] (const auto& haplotype_likelihoods) {
                               return haplotype_likelihoods.get()[i];
                           });
            
            result += Maths::log_sum_exp(tmp) - ln_ploidy_;
        }
        
        return result;
    }
}
} // namespace Octopus
