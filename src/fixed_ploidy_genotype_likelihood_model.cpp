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
    FixedPloidyGenotypeLikelihoodModel::FixedPloidyGenotypeLikelihoodModel(unsigned ploidy, const HaplotypeLikelihoodCache& haplotype_likelihoods)
    :
    haplotype_likelihoods_ {haplotype_likelihoods},
    ploidy_ {ploidy},
    ln_ploidy_ {std::log(ploidy)}
    {}
    
    // ln p(read | genotype) = ln sum {haplotype in genotype} p(read | haplotype) - ln ploidy
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
        if (genotype.is_homozygous()) {
            const auto& log_likelihoods = haplotype_likelihoods_.get().log_likelihoods(sample, genotype[0]);
            return std::accumulate(std::cbegin(log_likelihoods), std::cend(log_likelihoods), 0.0);
        }
        
        const auto& log_likelihoods1 = haplotype_likelihoods_.get().log_likelihoods(sample, genotype[0]);
        const auto& log_likelihoods2 = haplotype_likelihoods_.get().log_likelihoods(sample, genotype[1]);
        
        return std::inner_product(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1),
                                  std::cbegin(log_likelihoods2), 0.0, std::plus<double> {},
                                  [this] (const double a, const double b) -> double {
                                      return Maths::log_sum_exp(a, b) - ln_ploidy_;
                                  });
    }
    
    double FixedPloidyGenotypeLikelihoodModel::log_likelihood_triploid(const SampleIdType& sample,
                                                                        const Genotype<Haplotype>& genotype) const
    {
//        const auto& log_likelihoods1 = haplotype_likelihoods_.get().log_likelihoods(sample, genotype[0]);
//        const auto& log_likelihoods2 = haplotype_likelihoods_.get().log_likelihoods(sample, genotype[1]);
//        const auto& log_likelihoods3 = haplotype_likelihoods_.get().log_likelihoods(sample, genotype[1]);
        
        // TODO
        return 0;
//        return Maths::log_sum_exp(log_likelihood(read, genotype[0]),
//                                  log_likelihood(read, genotype[1]),
//                                  log_likelihood(read, genotype[2])) - ln_ploidy_;
    }
    
    double FixedPloidyGenotypeLikelihoodModel::log_likelihood_polyploid(const SampleIdType& sample,
                                                                         const Genotype<Haplotype>& genotype) const
    {
        // TODO
        return 0;
//        std::vector<double> log_haplotype_probabilities(ploidy_);
//        
//        std::transform(std::cbegin(genotype), std::cend(genotype),
//                       std::begin(log_haplotype_probabilities),
//                       [this, &read] (const auto& haplotype) {
//                           return log_likelihood(read, haplotype);
//                       });
//        
//        return Maths::log_sum_exp<double>(log_haplotype_probabilities) - ln_ploidy_;
    }
    
    // non-member methods
    
    namespace debug
    {
//        void print_read_genotype_liklihoods(const std::vector<Genotype<Haplotype>>& genotypes,
//                                            const ReadMap& reads,
//                                            const FixedPloidyGenotypeLikelihoodModel& read_model,
//                                            const size_t n)
//        {
//            auto m = std::min(n, genotypes.size());
//            
//            std::cout << "top " << n << " genotype likelihoods for each read in each sample" << std::endl;
//            
//            for (const auto& sample_reads : reads) {
//                std::cout << "Sample: " << sample_reads.first << ":" << std::endl;
//                for (const auto& read : sample_reads.second) {
//                    std::cout << "\tRead: " << read.get_region() << " " << read.get_cigar_string() << ":" << std::endl;
//                    
//                    std::vector<std::pair<Genotype<Haplotype>, double>> top {};
//                    top.reserve(genotypes.size());
//                    
//                    for (const auto& genotype : genotypes) {
//                        top.emplace_back(genotype, read_model.log_likelihood(read, genotype));
//                    }
//                    
//                    std::sort(std::begin(top), std::end(top),
//                              [] (const auto& lhs, const auto& rhs) {
//                                  return lhs.second > rhs.second;
//                              });
//                    
//                    for (unsigned i {}; i < m; ++i) {
//                        std::cout << "\t\t* ";
//                        print_variant_alleles(top[i].first);
//                        std::cout << " " << std::setprecision(10) << top[i].second << std::endl;
//                    }
//                }
//            }
//        }
    } // namespace debug
}
} // namespace Octopus
