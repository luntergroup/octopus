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
#include <iostream>

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
    
    double FixedPloidyGenotypeLikelihoodModel::log_probability(const AlignedRead& read, const Haplotype& haplotype) const
    {
        return haplotype_likelihoods_.get().log_probability(read, haplotype);
    }
    
    // ln p(read | genotype) = ln sum {haplotype in genotype} p(read | haplotype) - ln ploidy
    double FixedPloidyGenotypeLikelihoodModel::log_probability(const AlignedRead& read, const Genotype<Haplotype>& genotype) const
    {
        // These cases are just for optimisation
        switch (ploidy_) {
            case 1:
                return log_probability_haploid(read, genotype);
            case 2:
                return log_probability_diploid(read, genotype);
            case 3:
                return log_probability_triploid(read, genotype);
            default:
                return log_probability_polyploid(read, genotype);
        }
    }
    
    double FixedPloidyGenotypeLikelihoodModel::log_probability_haploid(const AlignedRead& read, const Genotype<Haplotype>& genotype) const
    {
        return log_probability(read, genotype.at(0));
    }
    
    double FixedPloidyGenotypeLikelihoodModel::log_probability_diploid(const AlignedRead& read, const Genotype<Haplotype>& genotype) const
    {
        return Maths::log_sum_exp(log_probability(read, genotype.at(0)),
                                  log_probability(read, genotype.at(1))) - ln_ploidy_;
    }
    
    double FixedPloidyGenotypeLikelihoodModel::log_probability_triploid(const AlignedRead& read, const Genotype<Haplotype>& genotype) const
    {
        return Maths::log_sum_exp(log_probability(read, genotype.at(0)),
                                  log_probability(read, genotype.at(1)),
                                  log_probability(read, genotype.at(2))) - ln_ploidy_;
    }
    
    double FixedPloidyGenotypeLikelihoodModel::log_probability_polyploid(const AlignedRead& read, const Genotype<Haplotype>& genotype) const
    {
        std::vector<double> log_haplotype_probabilities(ploidy_);
        
        std::transform(std::cbegin(genotype), std::cend(genotype), std::begin(log_haplotype_probabilities),
                       [this, &read] (const auto& haplotype) {
                           return log_probability(read, haplotype);
                       });
        
        return Maths::log_sum_exp<double>(log_haplotype_probabilities) - ln_ploidy_;
    }
    
    // non-member methods
    
    ProbabilityMatrix<Genotype<Haplotype>>
    compute_genotype_log_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                     const ReadMap& reads,
                                     const FixedPloidyGenotypeLikelihoodModel& read_model)
    {
        ProbabilityMatrix<Genotype<Haplotype>> result {};
        
        return result;
    }
    
    namespace debug
    {
        void print_genotype_log_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                            const FixedPloidyGenotypeLikelihoodModel& read_model,
                                            const size_t n)
        {
            // TODO
        }
        
        void print_read_genotype_liklihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                            const ReadMap& reads,
                                            const FixedPloidyGenotypeLikelihoodModel& read_model,
                                            const size_t n)
        {
            auto m = std::min(n, genotypes.size());
            
            std::cout << "top " << n << " genotype likelihoods for each read in each sample" << std::endl;
            
            for (const auto& sample_reads : reads) {
                std::cout << "Sample: " << sample_reads.first << ":" << std::endl;
                for (const auto& read : sample_reads.second) {
                    std::cout << "\tRead: " << read.get_region() << " " << read.get_cigar_string() << ":" << std::endl;
                    
                    std::vector<std::pair<Genotype<Haplotype>, double>> top {};
                    top.reserve(genotypes.size());
                    
                    for (const auto& genotype : genotypes) {
                        top.emplace_back(genotype, read_model.log_probability(read, genotype));
                    }
                    
                    std::sort(std::begin(top), std::end(top),
                              [] (const auto& lhs, const auto& rhs) {
                                  return lhs.second > rhs.second;
                              });
                    
                    for (unsigned i {}; i < m; ++i) {
                        std::cout << "\t\t* ";
                        print_variant_alleles(top[i].first);
                        std::cout << " " << std::setprecision(10) << top[i].second << std::endl;
                    }
                }
            }
        }
    } // namespace debug
}
} // namespace Octopus
