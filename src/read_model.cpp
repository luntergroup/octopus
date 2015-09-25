//
//  read_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "read_model.hpp"

#include <cmath>
#include <algorithm> // std::max, std::min

#include "pair_hmm.hpp"
#include "maths.hpp"

#include <iostream> // TEST

namespace Octopus
{

ReadModel::ReadModel(unsigned ploidy, bool can_cache_reads)
:
ploidy_ {ploidy},
can_cache_reads_ {can_cache_reads},
read_log_probability_cache_ {},
genotype_log_probability_cache_ {},
ln_ploidy_ {std::log(ploidy)}
{}

ReadModel::RealType ReadModel::log_probability(const AlignedRead& read, const Haplotype& haplotype,
                                               SampleIdType sample)
{
    if (is_read_in_cache(sample, read, haplotype)) {
        return get_log_probability_from_cache(sample, read, haplotype);
    }
    
    // TODO: make these members when pair_hmm is finalised
    
    RandomModel<RealType> lhs_random {};
    lhs_random.target_emission_probability = 0.25;
    lhs_random.query_emission_probability  = 0.25;
    
    MatchModel<RealType> match {};
    match.match_probability      = 0.25;
    match.gap_open_probability   = 0.015; // TODO: should be part of an error model
    match.gap_extend_probability = 0.020; // TODO: should be part of an error model
    
    RandomModel<RealType> rhs_random {};
    rhs_random.target_emission_probability = 0.25;
    rhs_random.query_emission_probability  = 0.25;
    
    // m.end_probability must satisfy:
    // m.end_probability <= 1 - 2 * m.gap_open_probability
    // m.end_probability <= 1 - m.gap_extend_probability
    
    auto max_match_end_prob = 1 - std::max(2 * match.gap_open_probability, match.gap_extend_probability);
    
    if (overlaps(read, haplotype)) {
        auto overlapped_region = get_overlapped(read, haplotype);
        auto covered_region    = get_encompassing(read, haplotype);
        
        if (begins_before(read, haplotype)) {
            lhs_random.target_end_probability = 0.99;
            lhs_random.query_end_probability  = 1.0 / (size(get_left_overhang(covered_region, overlapped_region)) + 1);
        } else {
            lhs_random.target_end_probability = 1.0 / (size(get_left_overhang(covered_region, overlapped_region)) + 1);
            lhs_random.query_end_probability  = 0.99;
        }
        
        match.end_probability = std::min(1.0 / (size(overlapped_region) + 1), max_match_end_prob);
        
        if (ends_before(read, haplotype)) {
            rhs_random.target_end_probability = 1.0 / (size(get_right_overhang(covered_region, overlapped_region)) + 1);
            rhs_random.query_end_probability  = 0.99;
        } else {
            rhs_random.target_end_probability = 0.99;
            rhs_random.query_end_probability  = 1.0 / (size(get_right_overhang(covered_region, overlapped_region)) + 1);
        }
    } else {
        lhs_random.target_end_probability = 1.0 / (size(haplotype) + 1);
        lhs_random.query_end_probability  = 1.0 / (size(read) + 1);
        
        match.end_probability = max_match_end_prob;
        
        rhs_random.target_end_probability = 0.99;
        rhs_random.query_end_probability  = 0.99;
    }
    
    auto joint_log_probability = nuc_log_viterbi_local<RealType>(haplotype.get_sequence(), read.get_sequence(),
                                                                 read.get_qualities(),
                                                                 lhs_random, match, rhs_random);
    
    auto conditional_log_probability = joint_log_probability -
                    haplotype.get_sequence().size() * std::log(lhs_random.target_emission_probability);
    
    add_read_to_cache(sample, read, haplotype, conditional_log_probability);
    
    return conditional_log_probability;
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

bool ReadModel::is_read_in_cache(SampleIdType sample, const AlignedRead& read,
                                 const Haplotype& haplotype) const noexcept
{
    if (!can_cache_reads_) return false;
    if (read_log_probability_cache_.count(sample) == 0) return false;
    if (read_log_probability_cache_.at(sample).count(read) == 0) return false;
    return read_log_probability_cache_.at(sample).at(read).count(haplotype) > 0;
}

ReadModel::RealType ReadModel::get_log_probability_from_cache(SampleIdType sample, const AlignedRead& read,
                                                              const Haplotype& haplotype) const
{
    return read_log_probability_cache_.at(sample).at(read).at(haplotype);
}

void ReadModel::add_read_to_cache(SampleIdType sample, const AlignedRead& read, const Haplotype& haplotype,
                                  RealType read_log_probability)
{
    if (can_cache_reads_) {
        read_log_probability_cache_[sample][read][haplotype] = read_log_probability;
    }
}

void ReadModel::clear_cache()
{
    read_log_probability_cache_.clear();
    genotype_log_probability_cache_.clear();
}

} // end namespace Octopus
