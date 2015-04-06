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
#include "aligned_read.h"

class ReadModel
{
public:
    using Reads = std::vector<AlignedRead>;
    
    ReadModel() = delete;
    explicit ReadModel(unsigned ploidy, bool can_cache_reads = true);
    ~ReadModel() = default;
    
    ReadModel(const ReadModel&)            = default;
    ReadModel& operator=(const ReadModel&) = default;
    ReadModel(ReadModel&&)                 = default;
    ReadModel& operator=(ReadModel&&)      = default;
    
    // ln p(read | haplotype)
    double log_probability(const AlignedRead& read, const Haplotype& haplotype, unsigned sample);
    
    // ln p(read | genotype)
    double log_probability(const AlignedRead& read, const Genotype& genotype, unsigned sample);
    
    // ln p(reads | genotype)
    double log_probability(const Reads& reads, const Genotype& genotype, unsigned sample);
    
private:
    unsigned ploidy_;
    bool can_cache_reads_;
    std::unordered_map<unsigned, std::unordered_map<AlignedRead,
        std::unordered_map<Haplotype, double>>> read_log_probability_cache_;
    std::unordered_map<unsigned, std::unordered_map<Genotype, double>> genotype_log_probability_cache_;
    
    const double ln_ploidy_;
    
    bool is_read_in_cache(unsigned sample, const AlignedRead& read, const Haplotype& haplotype) const noexcept;
    bool is_genotype_in_cache(unsigned sample, const Genotype& genotype) const noexcept;
    void add_read_to_cache(unsigned sample, const AlignedRead& read, const Haplotype& haplotype,
                           double read_log_probability);
    void add_genotype_to_cache(unsigned sample, const Genotype& genotype, double genotype_log_probability);
    
    // These are just for optimisation
    double log_probability_haploid(const AlignedRead& read, const Genotype& genotype,
                                   unsigned sample);
    double log_probability_diploid(const AlignedRead& read, const Genotype& genotype,
                                   unsigned sample);
    double log_probability_triploid(const AlignedRead& read, const Genotype& genotype,
                                    unsigned sample);
    double log_probability_polyploid(const AlignedRead& read, const Genotype& genotype,
                                     unsigned sample);
};

#endif /* defined(__Octopus__read_model__) */
