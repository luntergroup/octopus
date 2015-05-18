//
//  read_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "read_model.h"

#include <cmath>
#include <algorithm> // std::max, std::min

#include "pair_hmm.h"
#include "maths.h"

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
    
    RandomModel<RealType> r1 {};
    r1.background_probability = 0.25;
    
    MatchModel<RealType> m {};
    m.match_probability      = 1.0;
    m.gap_open_probability   = 0.017;
    m.gap_extend_probability = 0.025;
    
    RandomModel<RealType> r2 {};
    r2.background_probability = 0.25;
    
    // m.end_probability must satisfy:
    // m.end_probability <= 1 - 2 * m.gap_open_probability
    // m.end_probability <= 1 - m.gap_extend_probability
    auto m_end_max = 1 - std::max(2 * m.gap_open_probability, m.gap_extend_probability);
    
    if (overlaps(read, haplotype)) {
        auto overlapped_region = get_overlapped(read, haplotype);
        auto covered_region    = get_encompassing(read, haplotype);
        
        r1.end_probability = 1.0 / (size(get_left_overhang(covered_region, overlapped_region)) + 1);
        m.end_probability  = std::min(1.0 / (size(overlapped_region) + 1), m_end_max);
        r2.end_probability = 1.0 / (size(get_right_overhang(covered_region, overlapped_region)) + 1);
    } else {
        r1.end_probability = 1.0 / (std::max(size(read), size(haplotype)) + 1);
        m.end_probability  = m_end_max;
        r2.end_probability = 0.99;
    }
    
    auto result = nuc_log_viterbi_local<RealType>(haplotype.get_sequence(), read.get_sequence(),
                                                  read.get_qualities(), m, r1, r2);
    
    add_read_to_cache(sample, read, haplotype, result);
    
    return result;
}

// ln p(read | genotype) = ln sum {haplotype in genotype} p(read | haplotype) - ln ploidy
ReadModel::RealType ReadModel::log_probability(const AlignedRead& read, const Genotype<Haplotype>& genotype,
                                               SampleIdType sample)
{
    // These cases are just for optimisation; they are functionally equivalent
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
    //TODO
    std::vector<RealType> log_haplotype_probabilities {};
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
    
    auto log_sum_haplotype_probabilities = log_sum_exp<RealType>(log_haplotype_probabilities.cbegin(),
                                                               log_haplotype_probabilities.cend());
    
    return log_sum_haplotype_probabilities - ln_ploidy_;
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
