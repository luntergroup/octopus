//
//  haplotype_phaser.cpp
//  Octopus
//
//  Created by Daniel Cooke on 30/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "haplotype_phaser.h"

#include <algorithm>

#include "region_algorithms.h"
#include "allele.h"
#include "search_regions.h"
#include "haplotype.h"
#include "genotype.h"
#include "haplotype_prior_model.h"

HaplotypePhaser::HaplotypePhaser(ReferenceGenome& the_reference, VariationalBayesGenotypeModel& the_model,
                                 unsigned ploidy)
:
the_reference_ {the_reference},
ploidy_ {ploidy},
the_tree_ {the_reference},
the_model_ {the_model},
max_haplotypes_ {256},
max_region_density_ {8},
max_indicators_ {3},
max_model_update_iterations_ {3},
the_phased_genotypes_ {},
the_last_unphased_region_ {}
{}

HaplotypePhaser::PhasedGenotypePosteriors
HaplotypePhaser::get_phased_genotypes_posteriors(bool include_partially_phased_haplotypes)
{
    if (include_partially_phased_haplotypes) {
        remove_phased_region(the_last_unphased_region_);
    }
    
    return the_phased_genotypes_;
}

void HaplotypePhaser::phase()
{
    auto sub_region = the_last_unphased_region_;
    
    SamplesReads reads {};
    
    while (!the_candidates_.empty()) {
        auto candidate_range = overlap_range(the_candidates_.cbegin(), the_candidates_.cend(), sub_region);
        
        std::for_each(candidate_range.first, candidate_range.second, [this] (const Variant& candidate) {
            the_tree_.extend(candidate.get_reference_allele());
            the_tree_.extend(candidate.get_alternative_allele());
        });
        
        auto haplotypes = the_tree_.get_haplotypes(sub_region);
        
        unique_least_complex(haplotypes);
        
        auto genotypes = get_all_genotypes(haplotypes, ploidy_);
        
        auto haplotype_priors = get_haplotype_prior_probabilities<RealType>(haplotypes, candidate_range.first,
                                                                            candidate_range.second);
        
        Haplotype reference {the_reference_, sub_region};
        
        auto haplotype_pseudo_counts = get_prior_pseudo_counts(haplotype_priors, reference, 1.0);
        
        for (const auto& sample_reads : the_reads_) {
            reads.emplace(sample_reads.first, overlap_range(sample_reads.second.cbegin(),
                                                            sample_reads.second.cend(), sub_region));
        }
        
        auto genotype_posteriors = update_parameters(the_model_, genotypes, haplotype_pseudo_counts,
                                                     reads, max_model_update_iterations_);
        
        if (candidate_range.second != the_candidates_.cend()) {
            for (const auto& haplotype : haplotypes) {
                if (genotype_posteriors.second.at(haplotype) < haplotype_priors.at(haplotype) + 0.01) {
                    the_tree_.prune_all(haplotype);
                } else {
                    the_tree_.prune_unique(haplotype);
                }
            }
            
            return;
        }
        
        if (has_shared(the_reads_, sub_region, *candidate_range.second)) {
            for (const auto& haplotype : haplotypes) {
                if (genotype_posteriors.second.at(haplotype) < haplotype_priors.at(haplotype) + 0.01) {
                    the_tree_.prune_all(haplotype);
                } else {
                    the_tree_.prune_unique(haplotype);
                }
            }
            
            
            
            sub_region = next_sub_region(sub_region, the_reads_, the_candidates_, max_region_density_, max_indicators_);
        } else {
            the_phased_genotypes_.emplace_back(genotype_posteriors);
            
            the_candidates_.erase(the_candidates_.begin(), candidate_range.second);
            
            for (auto& sample_reads : the_reads_) {
                auto overlapped_reads = overlap_range(sample_reads.second.cbegin(),
                                                      sample_reads.second.cend(), *candidate_range.first);
                sample_reads.second.erase(sample_reads.second.begin(), overlapped_reads.first);
            }
            
            the_tree_.clear();
            the_model_.clear_cache();
            
            sub_region = next_sub_region(sub_region, the_reads_, the_candidates_, max_region_density_, max_indicators_);
        }
    }
}

void HaplotypePhaser::remove_phased_region(const GenomicRegion &a_region)
{
    
}
