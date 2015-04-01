//
//  read_model.h
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_model__
#define __Octopus__read_model__

#include <vector>
#include <unordered_map>

#include "haplotype.h"
#include "genotype.h"

class AlignedRead;

class ReadModel
{
public:
    using Reads = std::vector<AlignedRead>;
    
    ReadModel() = delete;
    explicit ReadModel(unsigned ploidy);
    ~ReadModel() = default;
    
    ReadModel(const ReadModel&)            = default;
    ReadModel& operator=(const ReadModel&) = default;
    ReadModel(ReadModel&&)                 = default;
    ReadModel& operator=(ReadModel&&)      = default;
    
    // ln p(read | haplotype)
    double log_probability(const AlignedRead& read, const Haplotype& haplotype) const;
    
    // ln p(reads | haplotype)
    double log_probability(const Reads& reads, const Haplotype& haplotype) const;
    
    // ln p(reads | genotype)
    double log_probability(const Reads& reads, const Genotype& genotype, unsigned sample);
    
private:
    unsigned ploidy_;
    std::unordered_map<unsigned, std::unordered_map<Haplotype, double>> haplotype_log_probability_cache_;
    
    const double ln_ploidy_;
    
    bool is_haplotype_in_cache(unsigned sample, const Haplotype& haplotype) const;
    
    // These are just for optimisation
    double log_probability_haploid(const Reads& reads, const Genotype& genotype,
                                   unsigned sample);
    double log_probability_diploid(const Reads& reads, const Genotype& genotype,
                                   unsigned sample);
    double log_probability_triploid(const Reads& reads, const Genotype& genotype,
                                    unsigned sample);
    double log_probability_polyploid(const Reads& reads, const Genotype& genotype,
                                     unsigned sample);
};

#endif /* defined(__Octopus__read_model__) */
