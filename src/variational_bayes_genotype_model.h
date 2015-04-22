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
    using RealType                       = double;
    using Haplotypes                     = std::vector<Haplotype>;
    using Genotypes                      = std::vector<Genotype>;
    using ReadIterator                   = ReadModel::ReadIterator;
    using HaplotypePseudoCounts          = std::unordered_map<Haplotype, RealType>;
    using SampleGenotypeResponsabilities = std::unordered_map<Genotype, RealType>;
    using GenotypeResponsabilities       = std::vector<SampleGenotypeResponsabilities>;
    
    VariationalBayesGenotypeModel() = delete;
    explicit VariationalBayesGenotypeModel(ReadModel& read_model, unsigned ploidy,
                                           RealType zero_epsilon=1e-20);
    ~VariationalBayesGenotypeModel() = default;
    
    VariationalBayesGenotypeModel(const VariationalBayesGenotypeModel&)            = default;
    VariationalBayesGenotypeModel& operator=(const VariationalBayesGenotypeModel&) = default;
    VariationalBayesGenotypeModel(VariationalBayesGenotypeModel&&)                 = default;
    VariationalBayesGenotypeModel& operator=(VariationalBayesGenotypeModel&&)      = default;
    
    RealType log_expected_genotype_probability(const Genotype& genotype,
                                               const HaplotypePseudoCounts& haplotype_pseudo_counts);
    
    RealType log_rho(const Genotype& genotype, const HaplotypePseudoCounts& haplotype_pseudo_counts,
                     ReadIterator first, ReadIterator last, unsigned sample);
    
    RealType genotype_responsability(const Genotype& genotype, ReadIterator first, ReadIterator last,
                                     const HaplotypePseudoCounts& haplotype_pseudo_counts,
                                     unsigned sample, const Genotypes& all_genotypes);
    
    SampleGenotypeResponsabilities genotype_responsabilities(const Genotypes& genotypes, ReadIterator first,
                                                             ReadIterator last,
                                                             const HaplotypePseudoCounts& haplotype_pseudo_counts,
                                                             unsigned sample);
    
    RealType expected_haplotype_count(const Haplotype& haplotype,
                                      const SampleGenotypeResponsabilities& sample_genotype_responsabilities);
    
    RealType posterior_haplotype_pseudo_count(const Haplotype& haplotype, RealType prior_pseudo_count,
                                              const GenotypeResponsabilities& genotype_responsabilities);
    
    RealType posterior_haplotype_probability(const Haplotype& haplotype,
                                             const HaplotypePseudoCounts& posterior_haplotype_pseudo_counts) const;
    
    RealType posterior_predictive_probability(const std::unordered_map<Haplotype, unsigned>& haplotype_counts,
                                              const HaplotypePseudoCounts& haplotype_pseudo_counts) const;
    
    RealType posterior_predictive_probability(const Genotype& genotype,
                                              const HaplotypePseudoCounts& haplotype_pseudo_counts) const;
    
    RealType posterior_probability_haplotype_in_samples(const Haplotype& haplotype,
                                                        const Genotypes& all_genotypes,
                                                        const HaplotypePseudoCounts& posterior_haplotype_pseudo_counts) const;
    
    RealType posterior_probability_haplotype_in_sample(const Haplotype& haplotype,
                                                       const Genotypes& all_genotypes,
                                                       const SampleGenotypeResponsabilities& genotype_responsabilities) const;
    
    RealType posterior_probability_allele_in_samples(const Allele& the_allele,
                                                     const Haplotypes& haplotypes,
                                                     const HaplotypePseudoCounts& posterior_haplotype_pseudo_counts) const;
    
    RealType posterior_probability_allele_in_sample(const Allele& the_allele,
                                                    const Haplotypes& haplotypes,
                                                    const SampleGenotypeResponsabilities& sample_genotype_responsabilities,
                                                    const Genotypes& genotypes) const;
    
private:
    unsigned ploidy_;
    ReadModel& read_model_;
    RealType zero_epsilon_;
    
    RealType sum(const HaplotypePseudoCounts& haplotype_pseudo_counts) const noexcept;
    
    // These are just for optimisation
    RealType log_expected_genotype_probability_haploid(const Genotype& genotype,
                                                       const HaplotypePseudoCounts& haplotype_pseudo_counts) const;
    RealType log_expected_genotype_probability_diploid(const Genotype& genotype,
                                                       const HaplotypePseudoCounts& haplotype_pseudo_counts) const;
    RealType log_expected_genotype_probability_triploid(const Genotype& genotype,
                                                        const HaplotypePseudoCounts& haplotype_pseudo_counts) const;
    RealType log_expected_genotype_probability_polyploid(const Genotype& genotype,
                                                         const HaplotypePseudoCounts& haplotype_pseudo_counts) const;
};

using GenotypePosteriors = std::pair<VariationalBayesGenotypeModel::GenotypeResponsabilities,
                                    VariationalBayesGenotypeModel::HaplotypePseudoCounts>;

using SamplesReads = std::vector<std::pair<VariationalBayesGenotypeModel::ReadIterator,
                                            VariationalBayesGenotypeModel::ReadIterator>>;

GenotypePosteriors update_parameters(VariationalBayesGenotypeModel& the_model,
                                     const VariationalBayesGenotypeModel::Genotypes& the_genotypes,
                                     const VariationalBayesGenotypeModel::HaplotypePseudoCounts& prior_haplotype_pseudocounts,
                                     const SamplesReads& the_reads, unsigned max_num_iterations);

#endif /* defined(__Octopus__variational_bayes_genotype_model.h__) */
