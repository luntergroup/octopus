//
//  empirical_variational_bayes_genotype_model.h
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__empirical_variational_bayes_genotype_model__
#define __Octopus__empirical_variational_bayes_genotype_model__

#include <vector>
#include <unordered_map>

#include "haplotype.h"
#include "genotype.h"
#include "read_model.h"

class AlignedRead;
class Variant;

class EmpiricalVariationalBayesGenotypeModel
{
public:
    using Haplotypes                     = std::vector<Haplotype>;
    using Genotypes                      = std::vector<Genotype>;
    using SampleReads                    = std::vector<AlignedRead>;
    using HaplotypePseudoCounts          = std::unordered_map<Haplotype, double>;
    using SampleGenotypeResponsabilities = std::unordered_map<Genotype, double>;
    using GenotypeResponsabilities       = std::vector<SampleGenotypeResponsabilities>;
    
    EmpiricalVariationalBayesGenotypeModel() = delete;
    explicit EmpiricalVariationalBayesGenotypeModel(ReadModel& read_model, unsigned ploidy);
    ~EmpiricalVariationalBayesGenotypeModel() = default;
    
    EmpiricalVariationalBayesGenotypeModel(const EmpiricalVariationalBayesGenotypeModel&)            = default;
    EmpiricalVariationalBayesGenotypeModel& operator=(const EmpiricalVariationalBayesGenotypeModel&) = default;
    EmpiricalVariationalBayesGenotypeModel(EmpiricalVariationalBayesGenotypeModel&&)                 = default;
    EmpiricalVariationalBayesGenotypeModel& operator=(EmpiricalVariationalBayesGenotypeModel&&)      = default;
    
    // E_pi [p(genotype | pi)]
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
                                            const HaplotypePseudoCounts& haplotype_pseudo_count) const;
    
    double allele_posterior_probability(const Variant& variant, const Haplotypes& haplotypes,
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

using GenotypePosteriors = std::pair<EmpiricalVariationalBayesGenotypeModel::GenotypeResponsabilities,
                                    EmpiricalVariationalBayesGenotypeModel::HaplotypePseudoCounts>;

using SamplesReads = std::vector<EmpiricalVariationalBayesGenotypeModel::SampleReads>;

GenotypePosteriors
update_parameters(EmpiricalVariationalBayesGenotypeModel& the_model,
                  const EmpiricalVariationalBayesGenotypeModel::Genotypes& the_genotypes,
                  const EmpiricalVariationalBayesGenotypeModel::HaplotypePseudoCounts& prior_haplotype_pseudocounts,
                  const SamplesReads& the_reads, unsigned max_num_iterations);

#endif /* defined(__Octopus__empirical_variational_bayes_genotype_model__) */
