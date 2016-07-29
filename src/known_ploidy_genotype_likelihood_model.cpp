//
//  known_ploidy_genotype_likelihood_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "known_ploidy_genotype_likelihood_model.hpp"

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <array>
#include <limits>
#include <cassert>

#include "maths.hpp"

#include <iostream> // TEST

namespace octopus { namespace model
{
static constexpr auto ln(const unsigned n)
{
    using T = double;
    
    constexpr std::array<T, 11> Ln {
        std::numeric_limits<T>::infinity(),
        0.0,
        0.693147180559945309417232121458176568075500134360255254120,
        1.098612288668109691395245236922525704647490557822749451734,
        1.386294361119890618834464242916353136151000268720510508241,
        1.609437912434100374600759333226187639525601354268517721912,
        1.791759469228055000812477358380702272722990692183004705855,
        1.945910149055313305105352743443179729637084729581861188459,
        2.079441541679835928251696364374529704226500403080765762362,
        2.197224577336219382790490473845051409294981115645498903469,
        2.302585092994045684017991454684364207601101488628772976033
    };
    return Ln[n];
}

KnownPloidyGenotypeLikelihoodModel::KnownPloidyGenotypeLikelihoodModel(unsigned ploidy,
                                                                       const HaplotypeLikelihoodCache& likelihoods)
:
likelihoods_ {likelihoods},
ploidy_ {ploidy}
{}

// ln p(read | genotype)  = ln sum {haplotype in genotype} p(read | haplotype) - ln ploidy
// ln p(reads | genotype) = sum {read in reads} ln p(read | genotype)
double KnownPloidyGenotypeLikelihoodModel::ln_likelihood(const Genotype<Haplotype>& genotype) const
{
    assert(likelihoods_.is_primed());
    
    // These cases are just for optimisation
    switch (ploidy_) {
        case 1:
            return ln_likelihood_haploid(genotype);
        case 2:
            return ln_likelihood_diploid(genotype);
        case 3:
            return ln_likelihood_triploid(genotype);
        case 4:
            return ln_likelihood_polyploid(genotype);
            //return log_likelihood_tetraploid(sample, genotype);
        default:
            return ln_likelihood_polyploid(genotype);
    }
}

double KnownPloidyGenotypeLikelihoodModel::ln_likelihood_haploid(const Genotype<Haplotype>& genotype) const
{
    const auto& log_likelihoods = likelihoods_[genotype[0]];
    return std::accumulate(std::cbegin(log_likelihoods), std::cend(log_likelihoods), 0.0);
}

double KnownPloidyGenotypeLikelihoodModel::ln_likelihood_diploid(const Genotype<Haplotype>& genotype) const
{
    const auto& log_likelihoods1 = likelihoods_[genotype[0]];
    
    if (genotype.is_homozygous()) {
        return std::accumulate(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1), 0.0);
    }
    
    const auto& log_likelihoods2 = likelihoods_[genotype[1]];
    
    return std::inner_product(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1),
                              std::cbegin(log_likelihoods2), 0.0, std::plus<> {},
                              [this] (const auto a, const auto b) -> double {
                                  return Maths::log_sum_exp(a, b) - ln(2);
                              });
}

double KnownPloidyGenotypeLikelihoodModel::ln_likelihood_triploid(const Genotype<Haplotype>& genotype) const
{
    using std::cbegin; using std::cend;
    
    const auto& log_likelihoods1 = likelihoods_[genotype[0]];
    
    if (genotype.is_homozygous()) {
        return std::accumulate(cbegin(log_likelihoods1), cend(log_likelihoods1), 0.0);
    }
    
    if (genotype.zygosity() == 3) {
        const auto& log_likelihoods2 = likelihoods_[genotype[1]];
        const auto& log_likelihoods3 = likelihoods_[genotype[2]];
        return Maths::inner_product(cbegin(log_likelihoods1), cend(log_likelihoods1),
                                    cbegin(log_likelihoods2), cbegin(log_likelihoods3),
                                    0.0, std::plus<> {},
                                    [this] (const auto a, const auto b, const auto c) -> double {
                                        return Maths::log_sum_exp(a, b, c) - ln(3);
                                    });
    }
    
    static const double ln2 {std::log(2)};
    
    if (genotype[0] != genotype[1]) {
        const auto& log_likelihoods2 = likelihoods_[genotype[1]];
        return std::inner_product(cbegin(log_likelihoods1), cend(log_likelihoods1),
                                  cbegin(log_likelihoods2), 0.0, std::plus<> {},
                                  [this] (const auto a, const auto b) -> double {
                                      return Maths::log_sum_exp(a, ln2 + b) - ln(3);
                                  });
    }
    
    const auto& log_likelihoods3 = likelihoods_[genotype[2]];
    
    return std::inner_product(cbegin(log_likelihoods1), cend(log_likelihoods1),
                              cbegin(log_likelihoods3), 0.0, std::plus<> {},
                              [this] (const auto a, const auto b) -> double {
                                  return Maths::log_sum_exp(ln2 + a, b) - ln(3);
                              });
}

double KnownPloidyGenotypeLikelihoodModel::ln_likelihood_tetraploid(const Genotype<Haplotype>& genotype) const
{
    const auto z = genotype.zygosity();
    
    const auto& log_likelihoods1 = likelihoods_[genotype[0]];
    
    if (z == 1) {
        return std::accumulate(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1), 0.0);
    }
    
    if (z == 4) {
        const auto& log_likelihoods2 = likelihoods_[genotype[1]];
        const auto& log_likelihoods3 = likelihoods_[genotype[2]];
        const auto& log_likelihoods4 = likelihoods_[genotype[3]];
        return Maths::inner_product(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1),
                                    std::cbegin(log_likelihoods2), std::cbegin(log_likelihoods3),
                                    std::cbegin(log_likelihoods4), 0.0, std::plus<> {},
                                    [this] (const auto a, const auto b, const auto c, const auto d) -> double {
                                        return Maths::log_sum_exp({a, b, c, d}) - ln(4);
                                    });
    }
    
    // TODO
    
    return 0;
}

double KnownPloidyGenotypeLikelihoodModel::ln_likelihood_polyploid(const Genotype<Haplotype>& genotype) const
{
    const auto z = genotype.zygosity();
    
    const auto& log_likelihoods1 = likelihoods_[genotype[0]];
    
    if (z == 1) {
        return std::accumulate(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1), 0.0);
    }
    
    if (z == 2) {
        static const double lnpm1 {std::log(ploidy_ - 1)};
        
        const auto unique_haplotypes = genotype.copy_unique_ref();
        
        const auto& log_likelihoods2 = likelihoods_[unique_haplotypes.back()];
        
        if (genotype.count(unique_haplotypes.front()) == 1) {
            return std::inner_product(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1),
                                      std::cbegin(log_likelihoods2), 0.0, std::plus<> {},
                                      [this] (const auto a, const auto b) -> double {
                                          return Maths::log_sum_exp(a, lnpm1 + b) - ln(ploidy_);
                                      });
        }
        
        return std::inner_product(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1),
                                  std::cbegin(log_likelihoods2), 0.0, std::plus<> {},
                                  [this] (const auto a, const auto b) -> double {
                                      return Maths::log_sum_exp(lnpm1 + a, b) - ln(ploidy_);
                                  });
    }
    
    std::vector<HaplotypeLikelihoodCache::LikelihoodVectorRef> ln_likelihoods {};
    ln_likelihoods.reserve(ploidy_);
    
    std::transform(std::cbegin(genotype), std::cend(genotype), std::back_inserter(ln_likelihoods),
                   [this] (const auto& haplotype)
                        -> const HaplotypeLikelihoodCache::LikelihoodVector& {
                       return likelihoods_[haplotype];
                   });
    
    std::vector<double> tmp(ploidy_);
    
    double result {0};
    
    const auto num_likelihoods = ln_likelihoods.front().get().size();
    
    for (std::size_t i {0}; i < num_likelihoods; ++i) {
        std::transform(std::cbegin(ln_likelihoods), std::cend(ln_likelihoods), std::begin(tmp),
                       [i] (const auto& likelihoods) {
                           return likelihoods.get()[i];
                       });
        
        result += Maths::log_sum_exp(tmp) - ln(ploidy_);
    }
    
    return result;
}
}
} // namespace octopus
