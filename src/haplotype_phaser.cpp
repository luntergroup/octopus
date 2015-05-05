//
//  haplotype_phaser.cpp
//  Octopus
//
//  Created by Daniel Cooke on 30/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "haplotype_phaser.h"

#include <algorithm> // std::for_each
#include <iterator>  // std::tie, std::distance
#include <cmath>     // std::log2

#include "allele.h"
#include "region_algorithms.h"
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
max_region_density_ {3},
max_model_update_iterations_ {3},
the_phased_genotypes_ {},
the_last_unphased_region_ {},
the_last_unphased_posteriors_ {}
{}

HaplotypePhaser::PhasedGenotypePosteriors
HaplotypePhaser::get_phased_genotypes_posteriors(bool include_partially_phased_haplotypes)
{
    if (include_partially_phased_haplotypes) {
        the_phased_genotypes_.emplace_back(the_last_unphased_posteriors_);
        remove_phased_region(the_candidates_.cbegin(), the_candidates_.cend());
    }
    
    return the_phased_genotypes_;
}

void HaplotypePhaser::phase()
{
    auto sub_region = the_last_unphased_region_;
    
    GenotypePosteriors genotype_posteriors;
    HaplotypePseudoCounts haplotype_prior_counts;
    HaplotypePseudoCounts haplotype_posterior_counts;
    
    unsigned max_candidates {max_region_density_};
    unsigned max_indicators {0};
    
    if (has_overlapped(the_candidates_.cbegin(), the_candidates_.cend(), sub_region)) {
        max_indicators = static_cast<unsigned>(count_overlapped(the_candidates_.cbegin(), the_candidates_.cend(),
                                                                sub_region));
        max_candidates = max_region_density_ + max_indicators;
        sub_region = next_sub_region(sub_region, the_reads_, the_candidates_, max_candidates, max_indicators);
    } else {
        sub_region = next_sub_region(sub_region, the_reads_, the_candidates_, max_candidates, max_indicators);
    }
    
    while (true) {
        cout << "looking in region " << sub_region << endl;
        
        auto candidate_range = overlap_range(the_candidates_.cbegin(), the_candidates_.cend(), sub_region);
        
        extend_haplotypes(candidate_range.first, candidate_range.second);
        
        auto haplotypes = get_haplotypes(sub_region);
        
        cout << "there are " << haplotypes.size() << " haplotypes" << endl;
        
        auto haplotype_prior_counts = get_haplotype_prior_counts(haplotypes, candidate_range.first,
                                                                 candidate_range.second, sub_region);
        
        std::tie(genotype_posteriors, haplotype_posterior_counts) =
            update_parameters(the_model_, get_all_genotypes(haplotypes, ploidy_), haplotype_prior_counts,
                              get_read_ranges(sub_region), max_model_update_iterations_);
        
        cout << "updated model" << endl;
        
        if (candidate_range.second == the_candidates_.cend()) {
            remove_unlikely_haplotypes(haplotypes, haplotype_prior_counts, haplotype_posterior_counts);
            the_last_unphased_region_ = sub_region;
            the_last_unphased_posteriors_ = {genotype_posteriors, haplotype_posterior_counts};
            return;
        }
        
        if (has_shared(the_reads_, sub_region, *candidate_range.second)) {
            remove_unlikely_haplotypes(haplotypes, haplotype_prior_counts, haplotype_posterior_counts);
            max_indicators = static_cast<unsigned>(std::distance(the_candidates_.cbegin(), candidate_range.second));
            unsigned k = std::log2(the_tree_.num_haplotypes());
            max_candidates = max_region_density_ + max_indicators - k;
            cout << "can phase with maximum " << max_indicators << " indicators and " << max_candidates << " candidates" << endl;
        } else {
            the_phased_genotypes_.emplace_back(std::make_pair(genotype_posteriors, haplotype_posterior_counts));
            remove_phased_region(candidate_range.first, candidate_range.second);
            max_indicators = 0;
            max_candidates = max_region_density_;
            cout << "could not phase" << endl;
        }
        
        cout << "there are now " << the_tree_.num_haplotypes() << " haplotypes" << endl;
        
        sub_region = next_sub_region(sub_region, the_reads_, the_candidates_, max_candidates, max_indicators);
    }
}

void HaplotypePhaser::extend_haplotypes(CandidateIterator first, CandidateIterator last)
{
    std::for_each(first, last, [this] (const Variant& candidate) {
        the_tree_.extend(candidate.get_reference_allele());
        the_tree_.extend(candidate.get_alternative_allele());
    });
}

HaplotypePhaser::Haplotypes HaplotypePhaser::get_haplotypes(const GenomicRegion& a_region)
{
    auto haplotypes = the_tree_.get_haplotypes(a_region);
    unique_least_complex(haplotypes);
    return haplotypes;
}

HaplotypePhaser::HaplotypePseudoCounts
HaplotypePhaser::get_haplotype_prior_counts(const Haplotypes& the_haplotypes,
                                            CandidateIterator first, CandidateIterator last,
                                            const GenomicRegion& the_region) const
{
    auto haplotype_priors = get_haplotype_prior_probabilities<RealType>(the_haplotypes, first, last);
    Haplotype reference {the_reference_, the_region};
    return get_haplotype_prior_pseudo_counts(haplotype_priors, reference, 1.0);
}

SamplesReads<HaplotypePhaser::ReadIterator> HaplotypePhaser::get_read_ranges(const GenomicRegion& the_region) const
{
    SamplesReads<ReadIterator> reads {};
    reads.reserve(the_reads_.size());
    
    for (const auto& sample_reads : the_reads_) {
        reads.emplace(sample_reads.first, overlap_range(sample_reads.second.cbegin(),
                                                        sample_reads.second.cend(), the_region));
    }
    
    return reads;
}

void HaplotypePhaser::remove_unlikely_haplotypes(const Haplotypes& the_haplotypes,
                                                 const HaplotypePseudoCounts& prior_counts,
                                                 const HaplotypePseudoCounts& posterior_counts)
{
    for (const auto& haplotype : the_haplotypes) {
        if (posterior_counts.at(haplotype) < prior_counts.at(haplotype) + 0.01) {
            the_tree_.prune_all(haplotype);
        } else {
            the_tree_.prune_unique(haplotype);
        }
    }
}

void HaplotypePhaser::remove_phased_region(CandidateIterator first, CandidateIterator last)
{
    the_candidates_.erase(the_candidates_.cbegin(), last);
    
    for (auto& sample_reads : the_reads_) {
        auto overlapped_reads = overlap_range(sample_reads.second.cbegin(),
                                              sample_reads.second.cend(), *first);
        sample_reads.second.erase(sample_reads.second.cbegin(), overlapped_reads.first);
    }
    
    the_tree_.clear();
    the_model_.clear_cache();
}
