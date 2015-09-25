//
//  read_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "read_model.hpp"

#include <cmath>     // std::log

#include "maths.hpp"

#include <iostream> // TEST

namespace Octopus
{

ReadModel::ReadModel(unsigned ploidy, bool can_cache_reads)
:
read_model_ {1000},
ploidy_ {ploidy},
can_cache_reads_ {can_cache_reads},
genotype_log_probability_cache_ {},
ln_ploidy_ {std::log(ploidy)}
{}

ReadModel::RealType ReadModel::log_probability(const AlignedRead& read, const Haplotype& haplotype,
                                               SampleIdType sample)
{
    return read_model_.log_probability(read, haplotype);
}

// ln p(read | genotype) = ln sum {haplotype in genotype} p(read | haplotype) - ln ploidy
ReadModel::RealType ReadModel::log_probability(const AlignedRead& read, const Genotype<Haplotype>& genotype,
                                               SampleIdType sample)
{
    // These cases are just for optimisation
    switch (ploidy_) {
        case 1:
            return log_probability_haploid(read, genotype, sample);
        case 2:
            return log_probability_diploid(read, genotype, sample);
        case 3:
            return log_probability_triploid(read, genotype, sample);
        default:
            return log_probability_polyploid(read, genotype, sample);
    }
}

void ReadModel::clear_cache()
{
    genotype_log_probability_cache_.clear();
}

ReadModel::RealType ReadModel::log_probability_haploid(const AlignedRead& read, const Genotype<Haplotype>& genotype,
                                                       SampleIdType sample)
{
    return log_probability(read, genotype.at(0), sample);
}

ReadModel::RealType ReadModel::log_probability_diploid(const AlignedRead& read, const Genotype<Haplotype>& genotype,
                                                       SampleIdType sample)
{
    return log_sum_exp(log_probability(read, genotype.at(0), sample),
                       log_probability(read, genotype.at(1), sample)) - ln_ploidy_;
}

ReadModel::RealType ReadModel::log_probability_triploid(const AlignedRead& read, const Genotype<Haplotype>& genotype,
                                                        SampleIdType sample)
{
    return log_sum_exp(log_probability(read, genotype.at(0), sample),
                       log_probability(read, genotype.at(1), sample),
                       log_probability(read, genotype.at(2), sample)) - ln_ploidy_;
}

ReadModel::RealType ReadModel::log_probability_polyploid(const AlignedRead& read, const Genotype<Haplotype>& genotype,
                                                         SampleIdType sample)
{
    std::vector<RealType> log_haplotype_probabilities(ploidy_);
    
    std::transform(std::cbegin(genotype), std::cend(genotype), log_haplotype_probabilities.begin(),
                   [this, &read, &sample] (const Haplotype& haplotype) {
                       return log_probability(read, haplotype, sample);
                   });
    
    return log_sum_exp<RealType>(log_haplotype_probabilities.cbegin(),
                                 log_haplotype_probabilities.cend()) - ln_ploidy_;
}
    
} // end namespace Octopus
