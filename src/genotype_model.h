//
//  genotype_model.h
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__genotype_model__
#define __Octopus__genotype_model__

#include <vector>
#include <cstddef>
#include <unordered_map>

#include "haplotype.h"
#include "genotype.h"
#include "read_model.h"

using std::size_t;

class AlignedRead;

class GenotypeModel
{
public:
    using Haplotypes  = std::vector<Haplotype>;
    using Genotypes   = std::vector<Genotype>;
    using SampleReads = std::vector<AlignedRead>;
    using HaplotypeLogProbabilities = std::unordered_map<Haplotype, double>;
    
    GenotypeModel() = delete;
    explicit GenotypeModel(ReadModel& read_model, unsigned ploidy);
    ~GenotypeModel() = default;
    
    GenotypeModel(const GenotypeModel&)            = default;
    GenotypeModel& operator=(const GenotypeModel&) = default;
    GenotypeModel(GenotypeModel&&)                 = default;
    GenotypeModel& operator=(GenotypeModel&&)      = default;
    
    // ln p(genotype | population_haplotype_log_probabilities)
    // Note this assumes Hardy-Weinberg equilibrium
    double log_probability(const Genotype& genotype,
                           const HaplotypeLogProbabilities& sample_haplotype_log_probabilities) const;
    
    // ln p(reads, genotype | population_haplotype_log_probabilities)
    double log_probability(const SampleReads& reads, const Genotype& genotype,
                           const HaplotypeLogProbabilities& sample_haplotype_log_probabilities,
                           unsigned sample);
    
    // ln p(reads | population_haplotype_log_probabilities)
    double log_probability(const SampleReads& reads,
                           const HaplotypeLogProbabilities& sample_haplotype_log_probabilities,
                           unsigned sample, const Genotypes& all_genotypes);
    
    // p(genotype | reads, population_haplotype_log_probabilities)
    double genotype_posterior(const Genotype& genotype, const SampleReads& reads,
                              const HaplotypeLogProbabilities& sample_haplotype_log_probabilities,
                              unsigned sample, const Genotypes& all_genotypes);
    
private:
    unsigned ploidy_;
    ReadModel& read_model_;
    
    // These are just for optimisation
    double log_probability_haploid(const Genotype& genotype,
                                   const HaplotypeLogProbabilities& sample_haplotype_log_probabilities) const;
    double log_probability_diploid(const Genotype& genotype,
                                   const HaplotypeLogProbabilities& sample_haplotype_log_probabilities) const;
    double log_probability_triploid(const Genotype& genotype,
                                    const HaplotypeLogProbabilities& sample_haplotype_log_probabilities) const;
    double log_probability_polyploid(const Genotype& genotype,
                                     const HaplotypeLogProbabilities& sample_haplotype_log_probabilities) const;
};

using Reads = std::vector<GenotypeModel::SampleReads>;
using HaplotypeProbabilities = std::unordered_map<Haplotype, double>;
using SampleGenotypeProbabilities = std::vector<std::unordered_map<Haplotype, double>>;

std::pair<HaplotypeProbabilities, SampleGenotypeProbabilities>
get_haplotype_probabilities(GenotypeModel the_model, GenotypeModel::Genotypes the_genotypes,
                            const Reads& the_reads);

void
update_haplotype_probabilities(GenotypeModel::Genotypes the_genotypes, 
                               GenotypeModel::HaplotypeLogProbabilities& haplotype_log_probabilities,
                               const Reads& the_reads, GenotypeModel the_model);

#endif /* defined(__Octopus__genotype_model__) */
