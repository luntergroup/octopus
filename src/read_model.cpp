//
//  read_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "read_model.h"

#include <cmath>

#include "aligned_read.h"
#include "pair_hmm.h"
#include "maths.h"

#include <iostream> //TEST

ReadModel::ReadModel(unsigned ploidy)
:
ploidy_ {ploidy},
haplotype_log_probability_cache_ {},
ln_ploidy_ {std::log(ploidy)}
{}

double ReadModel::log_probability(const AlignedRead& read, const Haplotype& haplotype) const
{
    RandomModel r {};
    r.background_probability = 0.25;
    r.end_probability        = 0.1;
    
    MatchModel m {};
    m.match_probability      = 1.0;
    m.gap_open_probability   = 0.0001;
    m.gap_extend_probability = 0.0001;
    m.end_probability        = 0.01;
    
    return nuc_log_viterbi_local<float>(haplotype.get_sequence(), read.get_sequence(),
                                        read.get_qualities(), m, r);
}

// ln p(reads | haplotype) = sum {read in reads} ln p(read | haplotype)
double ReadModel::log_probability(const Reads& reads, const Haplotype& haplotype) const
{
    double result {0}; // ln p(1) = 0
    
    for (const auto& read : reads) {
        result += log_probability(read, haplotype); // ln {p(a) * p(b)} = ln p(a) + ln p(b)
    }
    
    return result;
}

// ln p(reads | genotype) = ln sum {haplotype in genotype} p(reads | haplotype) - ln |genotype|
double ReadModel::log_probability(const Reads& reads, const Genotype& genotype,
                                  unsigned sample)
{
    // These cases are just for optimisation; they are functionally equivalent
    switch (ploidy_) {
        case 1:
            return log_probability_haploid(reads, genotype, sample);
        case 2:
            return log_probability_diploid(reads, genotype, sample);
        case 3:
            return log_probability_triploid(reads, genotype, sample);
        default:
            return log_probability_polyploid(reads, genotype, sample);
    }
}

double ReadModel::log_probability_haploid(const Reads& reads, const Genotype& genotype,
                                          unsigned sample)
{
    const auto& haplotype = genotype.at(0);
    double log_haplotype_probability;
    
    if (is_haplotype_in_cache(sample, haplotype)) {
        log_haplotype_probability = haplotype_log_probability_cache_.at(sample).at(haplotype);
    } else {
        log_haplotype_probability = log_probability(reads, haplotype);
        haplotype_log_probability_cache_[sample][haplotype] = log_haplotype_probability;
    }
    
    return log_haplotype_probability - ln_ploidy_;
}

double ReadModel::log_probability_diploid(const Reads& reads, const Genotype& genotype,
                                          unsigned sample)
{
    const auto& haplotype1 = genotype.at(0);
    const auto& haplotype2 = genotype.at(1);
    double log_haplotype1_probability;
    double log_haplotype2_probability;
    
    if (is_haplotype_in_cache(sample, haplotype1)) {
        log_haplotype1_probability = haplotype_log_probability_cache_.at(sample).at(haplotype1);
    } else {
        log_haplotype1_probability = log_probability(reads, haplotype1);
        haplotype_log_probability_cache_[sample][haplotype1] = log_haplotype1_probability;
    }
    
    if (is_haplotype_in_cache(sample, haplotype2)) {
        log_haplotype2_probability = haplotype_log_probability_cache_.at(sample).at(haplotype2);
    } else {
        log_haplotype2_probability = log_probability(reads, haplotype2);
        haplotype_log_probability_cache_[sample][haplotype2] = log_haplotype1_probability;
    }
    
    return log_sum_exp(log_haplotype1_probability, log_haplotype2_probability) - ln_ploidy_;
}

double ReadModel::log_probability_triploid(const Reads& reads, const Genotype& genotype,
                                           unsigned sample)
{
    const auto& haplotype1 = genotype.at(0);
    const auto& haplotype2 = genotype.at(1);
    const auto& haplotype3 = genotype.at(2);
    double log_haplotype1_probability;
    double log_haplotype2_probability;
    double log_haplotype3_probability;
    
    if (is_haplotype_in_cache(sample, haplotype1)) {
        log_haplotype1_probability = haplotype_log_probability_cache_.at(sample).at(haplotype1);
    } else {
        log_haplotype1_probability = log_probability(reads, haplotype1);
        haplotype_log_probability_cache_[sample][haplotype1] = log_haplotype1_probability;
    }
    
    if (is_haplotype_in_cache(sample, haplotype2)) {
        log_haplotype2_probability = haplotype_log_probability_cache_.at(sample).at(haplotype2);
    } else {
        log_haplotype2_probability = log_probability(reads, haplotype2);
        haplotype_log_probability_cache_[sample][haplotype2] = log_haplotype1_probability;
    }
    
    if (is_haplotype_in_cache(sample, haplotype3)) {
        log_haplotype3_probability = haplotype_log_probability_cache_.at(sample).at(haplotype3);
    } else {
        log_haplotype3_probability = log_probability(reads, haplotype3);
        haplotype_log_probability_cache_[sample][haplotype3] = log_haplotype1_probability;
    }
    
    return log_sum_exp(log_haplotype1_probability, log_haplotype2_probability,
                       log_haplotype3_probability) - ln_ploidy_;
}

double ReadModel::log_probability_polyploid(const Reads& reads, const Genotype& genotype,
                                            unsigned sample)
{
    std::vector<double> log_haplotype_probabilities {};
    log_haplotype_probabilities.reserve(ploidy_);
    
    //    for (const auto& haplotype : genotype) {
    //        if (is_haplotype_in_cache(sample, haplotype)) {
    //            log_haplotype_probabilities.push_back(haplotype_log_probability_cache_.at(sample).at(haplotype));
    //        } else {
    //            auto haplotype_log_probability = log_probability(reads, haplotype);
    //            haplotype_log_probability_cache_[sample][haplotype] = haplotype_log_probability;
    //            log_haplotype_probabilities.push_back(haplotype_log_probability);
    //        }
    //    }
    
    auto log_sum_haplotype_probabilities = log_sum_exp<double>(log_haplotype_probabilities.cbegin(),
                                                               log_haplotype_probabilities.cend());
    
    return log_sum_haplotype_probabilities - ln_ploidy_;
}

bool ReadModel::is_haplotype_in_cache(unsigned sample, const Haplotype& haplotype) const
{
    if (haplotype_log_probability_cache_.count(sample) == 0) return false;
    return haplotype_log_probability_cache_.at(sample).count(haplotype) > 0;
}
