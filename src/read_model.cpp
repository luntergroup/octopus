//
//  read_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "read_model.h"

#include <cmath>

#include "pair_hmm.h"
#include "maths.h"

#include <iostream> // TEST

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
    
    //TODO: make these members when pair_hmm is finalised
    RandomModel r {};
    r.background_probability = 0.25;
    r.end_probability        = 0.1;
    
    MatchModel m {};
    m.match_probability      = 1.0;
    m.gap_open_probability   = 0.017;
    m.gap_extend_probability = 0.001;
    m.end_probability        = 0.01;
    
    auto result = nuc_log_viterbi_local<float>(haplotype.get_sequence(), read.get_sequence(),
                                               read.get_qualities(), m, r);
    
    add_read_to_cache(sample, read, haplotype, result);
    
    return result;
}

// ln p(read | genotype) = ln sum {haplotype in genotype} p(read | haplotype) - ln ploidy
ReadModel::RealType ReadModel::log_probability(const AlignedRead& read, const Genotype& genotype,
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

// ln p(reads | genotype) = sum (read in reads} ln p(read | genotype)
ReadModel::RealType ReadModel::log_probability(ReadIterator first, ReadIterator last, const Genotype& genotype,
                                               SampleIdType sample)
{
    if (is_genotype_in_cache(sample, first, last, genotype)) {
        return get_log_probability_from_cache(sample, first, last, genotype);
    }
    
    RealType result {0};
    
    std::for_each(first, last, [this, &genotype, &sample, &result] (const auto& read) {
        result += log_probability(read, genotype, sample);
    });
    
    add_genotype_to_cache(sample, first, last, genotype, result);
    
    return result;
}

ReadModel::RealType ReadModel::log_probability_haploid(const AlignedRead& read, const Genotype& genotype,
                                                       SampleIdType sample)
{
    auto haplotype_log_probability = log_probability(read, genotype.at(0), sample);
    
    return haplotype_log_probability - ln_ploidy_;
}

ReadModel::RealType ReadModel::log_probability_diploid(const AlignedRead& read, const Genotype& genotype,
                                                       SampleIdType sample)
{
    auto haplotype1_log_probability = log_probability(read, genotype.at(0), sample);
    auto haplotype2_log_probability = log_probability(read, genotype.at(1), sample);
    
    return log_sum_exp(haplotype1_log_probability, haplotype2_log_probability) - ln_ploidy_;
}

ReadModel::RealType ReadModel::log_probability_triploid(const AlignedRead& read, const Genotype& genotype,
                                                        SampleIdType sample)
{
    auto haplotype1_log_probability = log_probability(read, genotype.at(0), sample);
    auto haplotype2_log_probability = log_probability(read, genotype.at(1), sample);
    auto haplotype3_log_probability = log_probability(read, genotype.at(2), sample);
    
    return log_sum_exp(haplotype1_log_probability, haplotype2_log_probability,
                       haplotype3_log_probability) - ln_ploidy_;
}

ReadModel::RealType ReadModel::log_probability_polyploid(const AlignedRead& read, const Genotype& genotype,
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

bool ReadModel::is_genotype_in_cache(SampleIdType sample, ReadIterator first, ReadIterator last,
                                     const Genotype& genotype) const noexcept
{
    if (genotype_log_probability_cache_.count(sample) == 0) return false;
    return genotype_log_probability_cache_.at(sample).count(std::make_tuple(genotype, *first, std::distance(first, last))) > 0;
}

void ReadModel::add_genotype_to_cache(SampleIdType sample, ReadIterator first, ReadIterator last,
                                      const Genotype& genotype, RealType genotype_log_probability)
{
    genotype_log_probability_cache_[sample][std::make_tuple(genotype, *first, std::distance(first, last))] = genotype_log_probability;
}

ReadModel::RealType ReadModel::get_log_probability_from_cache(SampleIdType sample, ReadIterator first,
                                                              ReadIterator last, const Genotype& genotype) const
{
    return genotype_log_probability_cache_.at(sample).at(std::make_tuple(genotype, *first, std::distance(first, last)));
}

void ReadModel::clear_cache()
{
    read_log_probability_cache_.clear();
    genotype_log_probability_cache_.clear();
}
