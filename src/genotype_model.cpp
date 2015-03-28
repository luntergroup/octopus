//
//  genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "genotype_model.h"

#include <cmath>

#include "maths.h"
#include "pair_hmm.h"

GenotypeModel::GenotypeModel(unsigned ploidy, size_t max_num_haplotypes)
:
ploidy_ {ploidy},
max_num_haplotypes_ {max_num_haplotypes},
sample_haplotype_log_probability_cache_ {},
ln_ploidy_ {std::log(ploidy)}
{}

std::vector<GenotypeModel::HaplotypePosteriors>
GenotypeModel::get_haplotype_probabilities(const Haplotypes& the_haplotypes,
                                           const SampleReads& the_reads)
{
    return std::vector<HaplotypePosteriors> {};
}

bool GenotypeModel::is_haplotype_in_cache(unsigned sample, const Haplotype& haplotype) const
{
    if (sample_haplotype_log_probability_cache_.count(sample) == 0) return false;
    return sample_haplotype_log_probability_cache_.at(sample).count(haplotype) > 0;
}

double GenotypeModel::log_probability(const AlignedRead& read, const Haplotype& haplotype) const
{
    RandomModel r {};
    r.background_probability = 0.25;
    r.end_probability = 0.6;
    
    MatchModel m {};
    m.match_probability = 1.0;
    m.gap_open_probability = 0.0001;
    m.gap_extend_probability = 0.0001;
    m.end_probability = 0.2;
    
    return nuc_log_viterbi_local<float>(haplotype.get_sequence(), read.get_sequence(),
                                        read.get_qualities(), m, r);
}

// ln p(reads | haplotype) = sum {read in reads} ln p(read | haplotype)
double GenotypeModel::log_probability(const Reads& reads, const Haplotype& haplotype) const
{
    double result {0}; // ln p(1) = 0
    
    for (const auto& read : reads) {
        result += log_probability(read, haplotype); // ln {p(a) * p(b)} = ln p(a) + ln p(b)
    }
    
    return result;
}

// ln p(reads | genotype) = ln sum {haplotype in genotype} p(reads | haplotype) - ln |genotype|
double GenotypeModel::log_probability(const Reads& reads, const Genotype& genotype,
                                      unsigned sample)
{
    // This seemingly bad polymorphic design is entirly for optimisation purposes; the
    // functionality in each case is entirly the same. Relax.
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

double GenotypeModel::log_probability(const Genotype& genotype,
                                      const HaplotypeLogProbabilities& population_haplotype_log_probabilities) const
{
    switch (ploidy_) {
        case 1:
            return log_probability_haploid(genotype, population_haplotype_log_probabilities);
        case 2:
            return log_probability_diploid(genotype, population_haplotype_log_probabilities);
        case 3:
            return log_probability_triploid(genotype, population_haplotype_log_probabilities);
        default:
            return log_probability_polyploid(genotype, population_haplotype_log_probabilities);
    }
}

// ln p(reads, genotype | php) = ln p(genotype | php) + ln p(reads | genotype)
// where php = population_haplotype_log_probabilities
double GenotypeModel::log_probability(const Reads& reads, const Genotype& genotype,
                                      const HaplotypeLogProbabilities& population_haplotype_log_probabilities,
                                      const unsigned sample)
{
    return log_probability(genotype, population_haplotype_log_probabilities) +
            log_probability(reads, genotype, sample);
}

double GenotypeModel::log_probability_haploid(const Reads& reads, const Genotype& genotype,
                                              unsigned sample)
{
    const auto& haplotype = genotype.front();
    double log_haplotype_probability;
    
    if (is_haplotype_in_cache(sample, haplotype)) {
        log_haplotype_probability = sample_haplotype_log_probability_cache_.at(sample).at(haplotype);
    } else {
        log_haplotype_probability = log_probability(reads, haplotype);
        sample_haplotype_log_probability_cache_[sample][haplotype] = log_haplotype_probability;
    }
    
    return log_haplotype_probability - ln_ploidy_;
}

double GenotypeModel::log_probability_diploid(const Reads& reads, const Genotype& genotype,
                                              unsigned sample)
{
    const auto& haplotype1 = genotype.front();
    const auto& haplotype2 = genotype.back();
    double log_haplotype1_probability;
    double log_haplotype2_probability;
    
    if (is_haplotype_in_cache(sample, haplotype1)) {
        log_haplotype1_probability = sample_haplotype_log_probability_cache_.at(sample).at(haplotype1);
    } else {
        log_haplotype1_probability = log_probability(reads, haplotype1);
        sample_haplotype_log_probability_cache_[sample][haplotype1] = log_haplotype1_probability;
    }
    
    if (is_haplotype_in_cache(sample, haplotype2)) {
        log_haplotype2_probability = sample_haplotype_log_probability_cache_.at(sample).at(haplotype2);
    } else {
        log_haplotype2_probability = log_probability(reads, haplotype2);
        sample_haplotype_log_probability_cache_[sample][haplotype2] = log_haplotype1_probability;
    }
    
    return log_sum_exp(log_haplotype1_probability, log_haplotype2_probability) - ln_ploidy_;
}

double GenotypeModel::log_probability_triploid(const Reads& reads, const Genotype& genotype,
                                               unsigned sample)
{
    const auto& haplotype1 = genotype.front();
    const auto& haplotype2 = genotype.at(1);
    const auto& haplotype3 = genotype.back();
    double log_haplotype1_probability;
    double log_haplotype2_probability;
    double log_haplotype3_probability;
    
    if (is_haplotype_in_cache(sample, haplotype1)) {
        log_haplotype1_probability = sample_haplotype_log_probability_cache_.at(sample).at(haplotype1);
    } else {
        log_haplotype1_probability = log_probability(reads, haplotype1);
        sample_haplotype_log_probability_cache_[sample][haplotype1] = log_haplotype1_probability;
    }
    
    if (is_haplotype_in_cache(sample, haplotype2)) {
        log_haplotype2_probability = sample_haplotype_log_probability_cache_.at(sample).at(haplotype2);
    } else {
        log_haplotype2_probability = log_probability(reads, haplotype2);
        sample_haplotype_log_probability_cache_[sample][haplotype2] = log_haplotype1_probability;
    }
    
    if (is_haplotype_in_cache(sample, haplotype3)) {
        log_haplotype3_probability = sample_haplotype_log_probability_cache_.at(sample).at(haplotype3);
    } else {
        log_haplotype3_probability = log_probability(reads, haplotype3);
        sample_haplotype_log_probability_cache_[sample][haplotype3] = log_haplotype1_probability;
    }
    
    return log_sum_exp(log_haplotype1_probability, log_haplotype2_probability,
                       log_haplotype3_probability) - ln_ploidy_;
}

double GenotypeModel::log_probability_polyploid(const Reads& reads, const Genotype& genotype,
                                                unsigned sample)
{
    std::vector<double> log_haplotype_probabilities {};
    log_haplotype_probabilities.reserve(ploidy_);
    
    for (const auto& haplotype : genotype) {
        if (is_haplotype_in_cache(sample, haplotype)) {
            log_haplotype_probabilities.push_back(sample_haplotype_log_probability_cache_.at(sample).at(haplotype));
        } else {
            auto haplotype_log_probability = log_probability(reads, haplotype);
            sample_haplotype_log_probability_cache_[sample][haplotype] = haplotype_log_probability;
            log_haplotype_probabilities.push_back(haplotype_log_probability);
        }
    }
    
    auto log_sum_haplotype_probabilities = log_sum_exp<double>(log_haplotype_probabilities.cbegin(),
                                                               log_haplotype_probabilities.cend());
    
    return log_sum_haplotype_probabilities - ln_ploidy_;
}

double GenotypeModel::log_probability_haploid(const Genotype& genotype,
                                              const HaplotypeLogProbabilities& population_haplotype_log_probabilities) const
{
    return population_haplotype_log_probabilities.at(genotype.front());
}

double GenotypeModel::log_probability_diploid(const Genotype& genotype,
                                              const HaplotypeLogProbabilities& population_haplotype_log_probabilities) const
{
    if (is_homozygous(genotype)) {
        return std::pow(population_haplotype_log_probabilities.at(genotype.front()), 2);
    } else {
        return 2 * population_haplotype_log_probabilities.at(genotype.front()) *
                    population_haplotype_log_probabilities.at(genotype.back());
    }
}

double GenotypeModel::log_probability_triploid(const Genotype& genotype,
                                               const HaplotypeLogProbabilities& population_haplotype_log_probabilities) const
{
    if (is_homozygous(genotype)) {
        return std::pow(population_haplotype_log_probabilities.at(genotype.front()), 3);
    } else {
        const auto& first_haplotype  = genotype.front();
        const auto& second_haplotype = get_first_alternate_haplotype(first_haplotype, genotype);
        if (num_occurences(first_haplotype, genotype) == 2) {
            return 3 * std::pow(population_haplotype_log_probabilities.at(first_haplotype), 2) *
                    population_haplotype_log_probabilities.at(second_haplotype);
        } else {
            if (num_occurences(second_haplotype, genotype) == 2) {
                return 3 * population_haplotype_log_probabilities.at(first_haplotype) *
                        std::pow(population_haplotype_log_probabilities.at(second_haplotype), 2);
            } else {
                return 6 * population_haplotype_log_probabilities.at(genotype.front()) *
                        population_haplotype_log_probabilities.at(genotype.at(1)) *
                        population_haplotype_log_probabilities.at(genotype.back());
            }
        }
    }
}

double GenotypeModel::log_probability_polyploid(const Genotype& genotype,
                                                const HaplotypeLogProbabilities& population_haplotype_log_probabilities) const
{
    return 0; // TODO: need multinomial
}
