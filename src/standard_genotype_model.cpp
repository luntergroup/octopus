//
//  standard_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "standard_genotype_model.h"

#include <cmath>

#include "aligned_read.h"
#include "maths.h"
#include "pair_hmm.h"

StandardGenotypeModel::StandardGenotypeModel(ReadModel& read_model, unsigned ploidy)
:
ploidy_ {ploidy},
read_model_ {read_model}
{}

// ln p(genotype | php) = sum {haplotype in genotype} (u(haplotype) * ln php[haplotype]) + ln C
// where u(haplotype) is the number of times the haplotype occurs in the genotypes
// and C is the constant multinomial coefficient described in documentation
double StandardGenotypeModel::log_probability(const Genotype& genotype,
                                      const HaplotypeLogProbabilities& sample_haplotype_log_probabilities) const
{
    // These cases are just for optimisation; they are functionally equivalent
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
double StandardGenotypeModel::log_probability(const SampleReads& reads, const Genotype& genotype,
                                      const HaplotypeLogProbabilities& sample_haplotype_log_probabilities,
                                      unsigned sample)
{
    return log_probability(genotype, sample_haplotype_log_probabilities) +
            read_model_.log_probability(reads.cbegin(), reads.cend(), genotype, sample);
}

double StandardGenotypeModel::log_probability(const SampleReads& reads,
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

double StandardGenotypeModel::genotype_posterior(const Genotype& genotype, const SampleReads& reads,
                                         const HaplotypeLogProbabilities& sample_haplotype_log_probabilities,
                                         unsigned sample, const Genotypes& all_genotypes)
{
    auto log_prior     = log_probability(genotype, sample_haplotype_log_probabilities);
    auto log_liklihood = read_model_.log_probability(reads.cbegin(), reads.cend(), genotype, sample);
    auto log_evidence  = log_probability(reads, sample_haplotype_log_probabilities, sample, all_genotypes);
    
    auto log_posterior = log_prior + log_liklihood - log_evidence; // Bayes theorem log form
    
    return std::exp(log_posterior);
}

double StandardGenotypeModel::log_probability_haploid(const Genotype& genotype,
                                              const HaplotypeLogProbabilities& sample_haplotype_log_probabilities) const
{
    static double ln_1 = std::log(1);
    return sample_haplotype_log_probabilities.at(genotype.at(0)) + ln_1;
}

double StandardGenotypeModel::log_probability_diploid(const Genotype& genotype,
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

double StandardGenotypeModel::log_probability_triploid(const Genotype& genotype,
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

double StandardGenotypeModel::log_probability_polyploid(const Genotype& genotype,
                                                const HaplotypeLogProbabilities& sample_haplotype_log_probabilities) const
{
    return 0; // TODO: need multinomial
}

std::pair<HaplotypeProbabilities, SampleGenotypeProbabilities>
get_haplotype_probabilities(StandardGenotypeModel the_model, StandardGenotypeModel::Genotypes the_genotypes,
                            const Reads& the_reads)
{
    return {HaplotypeProbabilities {}, SampleGenotypeProbabilities{}};
}

void update_haplotype_probabilities(StandardGenotypeModel::Genotypes the_genotypes,
                                    StandardGenotypeModel::HaplotypeLogProbabilities& haplotype_log_probabilities,
                                    const Reads& the_reads, StandardGenotypeModel the_model)
{
    const unsigned num_samples   = static_cast<unsigned>(the_reads.size());
    const unsigned num_genotypes = static_cast<unsigned>(the_genotypes.size());
    
    std::vector<std::vector<double>> genotype_responsabilities {};
    
    genotype_responsabilities.reserve(num_samples);
    
    for (unsigned n {0}; n < num_samples; ++n) {
    
        std::vector<double> sample_genotypes_responsabilities {};
        sample_genotypes_responsabilities.reserve(num_genotypes);
        
        for (const auto& genotype : the_genotypes) {
            
            sample_genotypes_responsabilities.push_back(the_model.genotype_posterior(genotype, the_reads.at(n),
                                                                             haplotype_log_probabilities, n, the_genotypes));
        }
        
        genotype_responsabilities.emplace_back(sample_genotypes_responsabilities);
    }
    
    
}
