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

GenotypeModel::GenotypeModel(unsigned ploidy)
:
ploidy_ {ploidy},
haplotype_log_probability_cache_ {},
ln_ploidy_ {std::log(ploidy)}
{}

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
double GenotypeModel::log_probability(const SampleReads& reads, const Haplotype& haplotype) const
{
    double result {0}; // ln p(1) = 0
    
    for (const auto& read : reads) {
        result += log_probability(read, haplotype); // ln {p(a) * p(b)} = ln p(a) + ln p(b)
    }
    
    return result;
}

// ln p(reads | genotype) = ln sum {haplotype in genotype} p(reads | haplotype) - ln |genotype|
double GenotypeModel::log_probability(const SampleReads& reads, const Genotype& genotype,
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

// ln p(genotype | php) = sum {haplotype in genotype} (u(haplotype) * ln php[haplotype]) + ln C
// where u(haplotype) is the number of times the haplotype occurs in the genotypes
// and C is the constant multinomial coefficient described in documentation
double GenotypeModel::log_probability(const Genotype& genotype,
                                      const HaplotypeLogProbabilities& sample_haplotype_log_probabilities) const
{
    switch (ploidy_) {
        case 1:
            return log_probability_haploid(genotype, sample_haplotype_log_probabilities);
        case 2:
            return log_probability_diploid(genotype, sample_haplotype_log_probabilities);
        case 3:
            return log_probability_triploid(genotype, sample_haplotype_log_probabilities);
        default:
            return log_probability_polyploid(genotype, sample_haplotype_log_probabilities);
    }
}

// ln p(reads, genotype | php) = ln p(genotype | php) + ln p(reads | genotype)
// where php = sample_haplotype_log_probabilities
double GenotypeModel::log_probability(const SampleReads& reads, const Genotype& genotype,
                                      const HaplotypeLogProbabilities& sample_haplotype_log_probabilities,
                                      unsigned sample)
{
    return log_probability(genotype, sample_haplotype_log_probabilities) +
            log_probability(reads, genotype, sample);
}

double GenotypeModel::log_probability(const SampleReads& reads,
                                      const HaplotypeLogProbabilities& sample_haplotype_log_probabilities,
                                      unsigned sample, const Genotypes& all_genotypes)
{
    std::vector<double> genotype_log_probabilities {};
    genotype_log_probabilities.reserve(all_genotypes.size());
    
    for (const auto& genotype : all_genotypes) {
        genotype_log_probabilities.push_back(log_probability(reads, genotype,
                                                             sample_haplotype_log_probabilities, sample));
    }
    
    return log_sum_exp<double>(genotype_log_probabilities.cbegin(), genotype_log_probabilities.cend());
}

double GenotypeModel::genotype_posterior(const Genotype& genotype, const SampleReads& reads,
                                         const HaplotypeLogProbabilities& sample_haplotype_log_probabilities,
                                         unsigned sample, const Genotypes& all_genotypes)
{
    auto log_prior     = log_probability(genotype, sample_haplotype_log_probabilities);
    auto log_liklihood = log_probability(reads, genotype, sample);
    auto log_evidence  = log_probability(reads, sample_haplotype_log_probabilities, sample, all_genotypes);
    
    auto log_posterior = log_prior + log_liklihood - log_evidence; // Bayes theorem log form
    
    return std::exp(log_posterior);
}

double GenotypeModel::log_probability_haploid(const SampleReads& reads, const Genotype& genotype,
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

double GenotypeModel::log_probability_diploid(const SampleReads& reads, const Genotype& genotype,
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

double GenotypeModel::log_probability_triploid(const SampleReads& reads, const Genotype& genotype,
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

double GenotypeModel::log_probability_polyploid(const SampleReads& reads, const Genotype& genotype,
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

double GenotypeModel::log_probability_haploid(const Genotype& genotype,
                                              const HaplotypeLogProbabilities& sample_haplotype_log_probabilities) const
{
    static double ln_1 = std::log(1);
    return sample_haplotype_log_probabilities.at(genotype.at(0)) + ln_1;
}

double GenotypeModel::log_probability_diploid(const Genotype& genotype,
                                              const HaplotypeLogProbabilities& sample_haplotype_log_probabilities) const
{
    static double ln_1 = std::log(1);
    static double ln_2 = std::log(2);
    if (genotype.is_homozygous()) {
        return 2 * sample_haplotype_log_probabilities.at(genotype.at(0)) + ln_1;
    } else {
        return sample_haplotype_log_probabilities.at(genotype.at(0)) +
                    sample_haplotype_log_probabilities.at(genotype.at(1)) + ln_2;
    }
}

double GenotypeModel::log_probability_triploid(const Genotype& genotype,
                                               const HaplotypeLogProbabilities& sample_haplotype_log_probabilities) const
{
//    static double ln_1 = std::log(1);
//    static double ln_3 = std::log(3);
//    static double ln_6 = std::log(6);
//    
//    if (genotype.is_homozygous()) {
//        return 3 * sample_haplotype_log_probabilities.at(genotype.at(0)) + ln_1;
//    } else {
//        const auto& first_haplotype  = genotype.at(0);
//        const auto& second_haplotype = genotype.get_first_alternate_haplotype(first_haplotype);
//        if (genotype.num_occurences(first_haplotype) == 2) {
//            return 2 * sample_haplotype_log_probabilities.at(first_haplotype) +
//                    sample_haplotype_log_probabilities.at(second_haplotype) + ln_3;
//        } else {
//            if (genotype.num_occurences(second_haplotype) == 2) {
//                return sample_haplotype_log_probabilities.at(first_haplotype) +
//                        2 * sample_haplotype_log_probabilities.at(second_haplotype) + ln_3;
//            } else {
//                return sample_haplotype_log_probabilities.at(genotype.at(0)) +
//                        sample_haplotype_log_probabilities.at(genotype.at(1)) +
//                        sample_haplotype_log_probabilities.at(genotype.at(2)) + ln_6;
//            }
//        }
//    }
    return 0;
}

double GenotypeModel::log_probability_polyploid(const Genotype& genotype,
                                                const HaplotypeLogProbabilities& sample_haplotype_log_probabilities) const
{
    return 0; // TODO: need multinomial
}

bool GenotypeModel::is_haplotype_in_cache(unsigned sample, const Haplotype& haplotype) const
{
    if (haplotype_log_probability_cache_.count(sample) == 0) return false;
    return haplotype_log_probability_cache_.at(sample).count(haplotype) > 0;
}

std::pair<HaplotypeProbabilities, SampleGenotypeProbabilities>
get_haplotype_probabilities(GenotypeModel the_model, GenotypeModel::Genotypes the_genotypes,
                            const Reads& the_reads)
{
    return {HaplotypeProbabilities {}, SampleGenotypeProbabilities{}};
}

void
update_haplotype_probabilities(GenotypeModel::Genotypes the_genotypes,
                               GenotypeModel::HaplotypeLogProbabilities& haplotype_log_probabilities,
                               const Reads& the_reads, GenotypeModel the_model)
{
    const unsigned num_samples   = static_cast<unsigned>(the_reads.size());
    const unsigned num_genotypes = static_cast<unsigned>(the_genotypes.size());
    
    std::vector<std::vector<double>> genotype_responsabilities {};
    
    genotype_responsabilities.reserve(num_samples);
    
    for (std::size_t n {0}; n < num_samples; ++n) {
    
        std::vector<double> sample_genotypes_responsabilities {};
        sample_genotypes_responsabilities.reserve(num_genotypes);
        
        for (const auto& genotype : the_genotypes) {
            
            sample_genotypes_responsabilities.push_back(the_model.genotype_posterior(genotype, the_reads.at(n),
                                                                             haplotype_log_probabilities, n, the_genotypes));
        }
        
        genotype_responsabilities.emplace_back(sample_genotypes_responsabilities);
    }
    
    
}
