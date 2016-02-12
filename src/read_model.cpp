//
//  read_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "read_model.hpp"

#include <vector>
#include <cmath>

#include "maths.hpp"

#include <iostream> // TEST

namespace Octopus
{

ReadModel::ReadModel(unsigned ploidy, const HaplotypeLikelihoodCache& haplotype_likelihoods)
:
haplotype_likelihoods_ {haplotype_likelihoods},
ploidy_ {ploidy},
ln_ploidy_ {std::log(ploidy)}
{}

double ReadModel::log_probability(const AlignedRead& read, const Haplotype& haplotype) const
{
    return haplotype_likelihoods_.get().log_probability(read, haplotype);
}

// ln p(read | genotype) = ln sum {haplotype in genotype} p(read | haplotype) - ln ploidy
double ReadModel::log_probability(const AlignedRead& read, const Genotype<Haplotype>& genotype) const
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

double ReadModel::log_probability_haploid(const AlignedRead& read, const Genotype<Haplotype>& genotype) const
{
    return log_probability(read, genotype.at(0));
}

double ReadModel::log_probability_diploid(const AlignedRead& read, const Genotype<Haplotype>& genotype) const
{
    return Maths::log_sum_exp(log_probability(read, genotype.at(0)),
                              log_probability(read, genotype.at(1))) - ln_ploidy_;
}

double ReadModel::log_probability_triploid(const AlignedRead& read, const Genotype<Haplotype>& genotype) const
{
    return Maths::log_sum_exp(log_probability(read, genotype.at(0)),
                              log_probability(read, genotype.at(1)),
                              log_probability(read, genotype.at(2))) - ln_ploidy_;
}

double ReadModel::log_probability_polyploid(const AlignedRead& read, const Genotype<Haplotype>& genotype) const
{
    std::vector<double> log_haplotype_probabilities(ploidy_);
    
    std::transform(std::cbegin(genotype), std::cend(genotype), std::begin(log_haplotype_probabilities),
                   [this, &read] (const auto& haplotype) {
                       return log_probability(read, haplotype);
                   });
    
    return Maths::log_sum_exp<double>(log_haplotype_probabilities) - ln_ploidy_;
}
    
} // namespace Octopus
