//
//  haplotype_phaser.cpp
//  Octopus
//
//  Created by Daniel Cooke on 30/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "haplotype_phaser.h"

#include <algorithm> // std::for_each
#include <iterator>  // std::distance
#include <cmath>     // std::log2
#include <utility>   // std::make_pair

#include "allele.h"
#include "region_algorithms.h"
#include "search_regions.h"
#include "haplotype.h"
#include "genotype.h"
#include "haplotype_prior_model.h"

HaplotypePhaser::HaplotypePhaser(ReferenceGenome& the_reference, VariationalBayesGenotypeModel& the_model,
                                 unsigned ploidy, unsigned max_haplotypes, unsigned max_model_update_iterations)
:
the_reference_ {the_reference},
ploidy_ {ploidy},
the_tree_ {the_reference},
the_model_ {the_model},
max_haplotypes_ {max_haplotypes},
max_region_density_ {static_cast<unsigned int>(std::log2(max_haplotypes_))},
max_model_update_iterations_ {max_model_update_iterations},
the_phased_regions_ {},
the_last_unphased_region_ {}
{}

HaplotypePhaser::PhasedRegions
HaplotypePhaser::get_phased_regions(bool include_partially_phased_regions)
{
    if (include_partially_phased_regions) {
        the_phased_regions_.emplace_back(std::move(the_last_unphased_region_));
        remove_phased_region(the_candidates_.cbegin(), the_candidates_.cend());
    }
    
    return the_phased_regions_;
}

void HaplotypePhaser::phase()
{
    HaplotypePseudoCounts haplotype_prior_counts;
    
    BayesianGenotypeModel::Latents<SampleIdType, RealType> latent_posteriors;
    
    GenomicRegion sub_region {the_last_unphased_region_.the_region};
    unsigned max_candidates {max_region_density_};
    unsigned max_indicators {0};
    std::pair<CandidateIterator, CandidateIterator> candidate_range;
    CandidateIterator previous_region_last_it {the_candidates_.cbegin()};
    
    if (has_overlapped(the_candidates_.cbegin(), the_candidates_.cend(), sub_region)) {
        max_indicators = static_cast<unsigned>(count_overlapped(the_candidates_.cbegin(), the_candidates_.cend(),
                                                                sub_region));
        max_candidates = max_region_density_ + max_indicators - std::log2(the_tree_.num_haplotypes());
        previous_region_last_it = overlap_range(the_candidates_.cbegin(), the_candidates_.cend(), sub_region).second;
    }
    
    sub_region = next_sub_region(sub_region, the_reads_, the_candidates_, max_candidates, max_indicators);
    
    while (true) {
        cout << "looking in region " << sub_region << endl;
        
        candidate_range = overlap_range(the_candidates_.cbegin(), the_candidates_.cend(), sub_region);
        
        extend_haplotypes(previous_region_last_it, candidate_range.second);
        
        auto haplotypes = get_haplotypes(sub_region);
        
        cout << "there are " << haplotypes.size() << " haplotypes" << endl;
        
        auto haplotype_prior_counts = get_haplotype_prior_counts(haplotypes, candidate_range.first,
                                                                 candidate_range.second, sub_region);
        
        auto genotypes = get_all_genotypes(haplotypes, ploidy_);
        
        latent_posteriors = BayesianGenotypeModel::update_latents(the_model_, genotypes, haplotype_prior_counts,
                                                                  get_read_ranges(sub_region), max_model_update_iterations_);
        
        cout << "updated model" << endl;
        
        if (candidate_range.second == the_candidates_.cend()) {
            remove_unlikely_haplotypes(haplotypes, haplotype_prior_counts, latent_posteriors.haplotype_pseudo_counts);
            the_last_unphased_region_ = PhasedRegion {sub_region, std::move(haplotypes),
                                                std::move(genotypes), std::move(latent_posteriors)};
            return;
        } else if (has_shared(the_reads_, sub_region, *candidate_range.second)) {
            remove_unlikely_haplotypes(haplotypes, haplotype_prior_counts, latent_posteriors.haplotype_pseudo_counts);
            max_indicators = static_cast<unsigned>(std::distance(the_candidates_.cbegin(), candidate_range.second));
            max_candidates = max_region_density_ + max_indicators - std::log2(the_tree_.num_haplotypes());
            previous_region_last_it = candidate_range.second;
            cout << "can phase with maximum " << max_indicators << " indicators and " << max_candidates << " candidates" << endl;
        } else {
            the_phased_regions_.emplace_back(sub_region, std::move(haplotypes), std::move(genotypes),
                                             std::move(latent_posteriors));
            remove_phased_region(candidate_range.first, candidate_range.second);
            max_indicators = 0;
            max_candidates = max_region_density_;
            previous_region_last_it = the_candidates_.cbegin();
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
    return BayesianGenotypeModel::get_haplotype_prior_pseudo_counts(haplotype_priors, reference, 1.0);
}

HaplotypePhaser::ReadRanges HaplotypePhaser::get_read_ranges(const GenomicRegion& the_region) const
{
    ReadRanges reads {};
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
        //haplotype.print_explicit_alleles();
        //cout << prior_counts.at(haplotype) << " " << posterior_counts.at(haplotype) << " " << (posterior_counts.at(haplotype) / prior_counts.at(haplotype)) << endl;
        if (posterior_counts.at(haplotype) / prior_counts.at(haplotype) < 1.001) {
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
