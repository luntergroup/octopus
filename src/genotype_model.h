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

using std::size_t;

class AlignedRead;

class GenotypeModel
{
public:
    using SampleProbabilities = std::vector<double>;
    struct HaplotypePosteriors
    {
        double population_probability;
        SampleProbabilities sample_probabilities;
    };
    using Haplotypes  = std::vector<Haplotype>;
    using SampleReads = std::vector<std::vector<AlignedRead>>;
    
    GenotypeModel() = delete;
    explicit GenotypeModel(unsigned ploidy, size_t max_num_haplotypes);
    ~GenotypeModel() = default;
    
    GenotypeModel(const GenotypeModel&)            = default;
    GenotypeModel& operator=(const GenotypeModel&) = default;
    GenotypeModel(GenotypeModel&&)                 = default;
    GenotypeModel& operator=(GenotypeModel&&)      = default;
    
    std::vector<HaplotypePosteriors> get_haplotype_probabilities(const Haplotypes& the_haplotypes,
                                                                 const SampleReads& the_reads);
private:
    using Reads                     = std::vector<AlignedRead>;
    using HaplotypeLogProbabilities = std::unordered_map<Haplotype, double>;
    
    unsigned ploidy_;
    size_t max_num_haplotypes_;
    std::unordered_map<unsigned, std::unordered_map<Haplotype, double>> sample_haplotype_log_probability_cache_;
    
    const double ln_ploidy_;
    
    bool is_haplotype_in_cache(unsigned sample, const Haplotype& haplotype) const;
    
    // ln p(read | haplotype)
    double log_probability(const AlignedRead& read, const Haplotype& haplotype) const;
    
    // ln p(reads | haplotype)
    double log_probability(const Reads& reads, const Haplotype& haplotype) const;
    
    // ln p(reads | genotype)
    double log_probability(const Reads& reads, const Genotype& genotype, unsigned sample);
    
    // ln p(genotype | population_haplotype_log_probabilities)
    // Note this assumes Hardy-Weinberg equilibrium
    double log_probability(const Genotype& genotype,
                           const HaplotypeLogProbabilities& population_haplotype_log_probabilities) const;
    
    // ln p(reads, genotype | population_haplotype_log_probabilities)
    double log_probability(const Reads& reads, const Genotype& genotype,
                           const HaplotypeLogProbabilities& population_haplotype_log_probabilities,
                           const unsigned sample);
    
    // These are just for optimisation
    double log_probability_haploid(const Reads& reads, const Genotype& genotype, unsigned sample);
    double log_probability_diploid(const Reads& reads, const Genotype& genotype, unsigned sample);
    double log_probability_triploid(const Reads& reads, const Genotype& genotype, unsigned sample);
    double log_probability_polyploid(const Reads& reads, const Genotype& genotype, unsigned sample);
    
    double log_probability_haploid(const Genotype& genotype,
                                   const HaplotypeLogProbabilities& population_haplotype_log_probabilities) const;
    double log_probability_diploid(const Genotype& genotype,
                                   const HaplotypeLogProbabilities& population_haplotype_log_probabilities) const;
    double log_probability_triploid(const Genotype& genotype,
                                    const HaplotypeLogProbabilities& population_haplotype_log_probabilities) const;
    double log_probability_polyploid(const Genotype& genotype,
                                     const HaplotypeLogProbabilities& population_haplotype_log_probabilities) const;
};

#endif /* defined(__Octopus__genotype_model__) */
