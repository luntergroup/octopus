//
//  variational_bayes_genotype_model.h.h
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__variational_bayes_genotype_model__
#define __Octopus__variational_bayes_genotype_model__

#include <vector>
#include <unordered_map>

#include "haplotype.h"
#include "genotype.h"
#include "read_model.h"

class AlignedRead;

class VariationalBayesGenotypeModel
{
public:
    using Haplotypes                     = std::vector<Haplotype>;
    using Genotypes                      = std::vector<Genotype>;
    using SampleReads                    = std::vector<AlignedRead>;
    using HaplotypePseudoCounts          = std::unordered_map<Haplotype, double>;
    using SampleGenotypeResponsabilities = std::unordered_map<Genotype, double>;
    using GenotypeResponsabilities       = std::vector<SampleGenotypeResponsabilities>;
    
    VariationalBayesGenotypeModel() = delete;
    explicit VariationalBayesGenotypeModel(ReadModel& read_model, unsigned ploidy);
    ~VariationalBayesGenotypeModel() = default;
    
    VariationalBayesGenotypeModel(const VariationalBayesGenotypeModel&)            = default;
    VariationalBayesGenotypeModel& operator=(const VariationalBayesGenotypeModel&) = default;
    VariationalBayesGenotypeModel(VariationalBayesGenotypeModel&&)                 = default;
    VariationalBayesGenotypeModel& operator=(VariationalBayesGenotypeModel&&)      = default;
    
    double log_expected_genotype_probability(const Genotype& genotype,
                                             const HaplotypePseudoCounts& haplotype_pseudo_counts);
    
    double log_rho(const Genotype& genotype, const HaplotypePseudoCounts& haplotype_pseudo_counts,
                   const SampleReads& reads, unsigned sample);
    
    double genotype_responsability(const Genotype& genotype, const SampleReads& reads,
                                   const HaplotypePseudoCounts& haplotype_pseudo_counts,
                                   unsigned sample, const Genotypes& all_genotypes);
    
    double expected_haplotype_count(const Haplotype& haplotype,
                                    const SampleGenotypeResponsabilities& sample_genotype_responsabilities);
    
    double posterior_haplotype_pseudo_count(const Haplotype& haplotype, double prior_pseudo_count,
                                            const GenotypeResponsabilities& genotype_responsabilities);
    
    double posterior_predictive_probability(const std::unordered_map<Haplotype, unsigned>& haplotype_counts,
                                            const HaplotypePseudoCounts& haplotype_pseudo_counts) const;
    
    double posterior_haplotype_probability(const Haplotype& haplotype,
                                           const HaplotypePseudoCounts& posterior_haplotype_pseudo_counts) const;
    
    double posterior_haplotype_probability(const Haplotype& haplotype,
                                           const SampleGenotypeResponsabilities& genotype_responsabilities) const;
    
    double allele_posterior_probability(const GenomicRegion& the_allele_region,
                                        const Haplotype::SequenceType& the_allele_sequence,
                                        const Haplotypes& haplotypes,
                                        const HaplotypePseudoCounts& posterior_haplotype_pseudo_counts) const;
    
    double allele_posterior_probability(const GenomicRegion& the_allele_region,
                                        const Haplotype::SequenceType& the_allele_sequence,
                                        const Haplotypes& haplotypes,
                                        const SampleGenotypeResponsabilities& sample_genotype_responsabilities,
                                        const Genotypes& genotypes) const;
    
private:
    unsigned ploidy_;
    ReadModel& read_model_;
    
    unsigned pseudo_count_sum(const HaplotypePseudoCounts& haplotype_pseudo_counts) const;
    
    // These are just for optimisation
    double log_expected_genotype_probability_haploid(const Genotype& genotype,
                                                     const HaplotypePseudoCounts& haplotype_pseudo_counts) const;
    double log_expected_genotype_probability_diploid(const Genotype& genotype,
                                                     const HaplotypePseudoCounts& haplotype_pseudo_counts) const;
    double log_expected_genotype_probability_triploid(const Genotype& genotype,
                                                      const HaplotypePseudoCounts& haplotype_pseudo_counts) const;
    double log_expected_genotype_probability_polyploid(const Genotype& genotype,
                                                       const HaplotypePseudoCounts& haplotype_pseudo_counts) const;
};

using GenotypePosteriors = std::pair<VariationalBayesGenotypeModel::GenotypeResponsabilities,
                                    VariationalBayesGenotypeModel::HaplotypePseudoCounts>;

using SamplesReads = std::vector<VariationalBayesGenotypeModel::SampleReads>;

GenotypePosteriors
update_parameters(VariationalBayesGenotypeModel& the_model,
                  const VariationalBayesGenotypeModel::Genotypes& the_genotypes,
                  const VariationalBayesGenotypeModel::HaplotypePseudoCounts& prior_haplotype_pseudocounts,
                  const SamplesReads& the_reads, unsigned max_num_iterations);

#endif /* defined(__Octopus__variational_bayes_genotype_model.h__) */
