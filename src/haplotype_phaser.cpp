//
//  haplotype_phaser.cpp
//  Octopus
//
//  Created by Daniel Cooke on 30/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "haplotype_phaser.hpp"

#include <algorithm> // std::for_each
#include <iterator>  // std::distance
#include <cmath>     // std::log2
#include <utility>   // std::make_pair

#include "allele.hpp"
#include "mappable_algorithms.hpp"
#include "search_regions.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "haplotype_prior_model.hpp"

namespace Octopus
{
//HaplotypePhaser::HaplotypePhaser(ReferenceGenome& the_reference, VariationalBayesGenotypeModel& the_model,
//                                 unsigned ploidy, unsigned max_haplotypes, unsigned max_model_update_iterations,
//                                 RealType reference_bias)
//:
//the_reference_ {the_reference},
//ploidy_ {ploidy},
//the_tree_ {the_reference},
//the_model_ {the_model},
//reference_bias_ {reference_bias},
//max_haplotypes_ {max_haplotypes},
//max_region_density_ {static_cast<unsigned int>(std::log2(max_haplotypes_))},
//max_model_update_iterations_ {max_model_update_iterations},
//the_phased_regions_ {},
//the_last_unphased_region_ {}
//{}
//
//HaplotypePhaser::PhasedRegions
//HaplotypePhaser::get_phased_regions(SteamingStatus status)
//{
//    if (status == SteamingStatus::Finished) {
//        the_phased_regions_.emplace_back(std::move(the_last_unphased_region_));
//        clear_phased_region(the_candidates_.cbegin(), the_candidates_.cend());
//    }
//    
//    return the_phased_regions_;
//}
//
//void HaplotypePhaser::phase_current_data()
//{
//    HaplotypePseudoCounts haplotype_prior_counts {};
//    BayesianGenotypeModel::Latents<SampleIdType, RealType> latent_posteriors {};
//    GenomicRegion sub_region {the_last_unphased_region_.the_region};
//    unsigned max_candidates {max_region_density_};
//    unsigned max_indicators {0};
//    CandidateIterator previous_region_last_candidate_it {the_candidates_.cbegin()};
//    
//    auto num_candidates_unphased = count_overlapped(the_candidates_.cbegin(), the_candidates_.cend(),
//                                                    sub_region);
//    
//    if (num_candidates_unphased > 0) {
//        max_indicators = static_cast<unsigned>(num_candidates_unphased);
//        max_candidates = max_region_density_ + max_indicators - std::log2(the_tree_.num_haplotypes());
//        previous_region_last_candidate_it = overlap_range(the_candidates_.cbegin(), the_candidates_.cend(),
//                                                          sub_region).end().base();
//    }
//    
//    sub_region = advance_region(sub_region, the_reads_, the_candidates_, max_candidates, max_indicators,
//                                IndicatorLimit::NoLimit, ExtensionLimit::NoLimit);
//    
//    while (true) {
////        cout << "looking in region " << sub_region << endl;
////        auto candidate_range = contained_range(the_candidates_.cbegin(), the_candidates_.cend(), sub_region);
////        
////        extend_tree(previous_region_last_candidate_it, candidate_range.end());
////        
////        auto haplotypes = get_unique_haplotypes_from_tree(sub_region);
////        auto haplotype_prior_counts = get_haplotype_prior_counts(haplotypes, candidate_range.begin(),
////                                                                 candidate_range.end(), sub_region);
////        auto genotypes = get_all_genotypes(haplotypes, ploidy_);
////        
////        cout << "proceeding with " << haplotypes.size() << " haplotypes" << endl;
////        
////        latent_posteriors = BayesianGenotypeModel::update_latents(the_model_, genotypes,
////                                                                  haplotype_prior_counts,
////                                                                  get_read_iterator_ranges(sub_region),
////                                                                  max_model_update_iterations_);
////        
////        cout << "computed latents" << endl;
////        
////        if (candidate_range.end() == the_candidates_.cend()) {
////            remove_unlikely_haplotypes_from_tree(haplotypes, haplotype_prior_counts,
////                                                 latent_posteriors.haplotype_pseudo_counts);
////            the_last_unphased_region_ = PhasedRegion {sub_region, std::move(haplotypes),
////                                                std::move(genotypes), std::move(latent_posteriors)};
////            return;
////        } else if (is_potentially_phasable(sub_region, *candidate_range.end())) {
////            remove_unlikely_haplotypes_from_tree(haplotypes, haplotype_prior_counts,
////                                                 latent_posteriors.haplotype_pseudo_counts);
////            
////            max_indicators = static_cast<unsigned>(std::distance(the_candidates_.cbegin(), candidate_range.end()));
////            max_candidates = max_region_density_ + max_indicators - std::log2(the_tree_.num_haplotypes());
////            previous_region_last_candidate_it = candidate_range.end();
////            cout << "can phase with maximum " << max_indicators << " indicators and " << max_candidates << " candidates" << endl;
////        } else {
////            the_phased_regions_.emplace_back(sub_region, std::move(haplotypes), std::move(genotypes),
////                                             std::move(latent_posteriors));
////            clear_phased_region(candidate_range.begin(), candidate_range.end());
////            
////            max_indicators = 0;
////            max_candidates = max_region_density_;
////            previous_region_last_candidate_it = the_candidates_.cbegin();
////            cout << "could not phase" << endl;
////        }
////        
////        cout << "there are now " << the_tree_.num_haplotypes() << " haplotypes" << endl;
////        
////        sub_region = next_sub_region(sub_region, the_reads_, the_candidates_, max_candidates, max_indicators,
////                                     IndicatorLimit::NoLimit, ExtensionLimit::NoLimit);
//    }
//}
//
//void HaplotypePhaser::extend_tree(CandidateIterator first, CandidateIterator last)
//{
//    std::for_each(first, last, [this] (const Variant& candidate) {
//        the_tree_.extend(candidate.get_reference_allele());
//        the_tree_.extend(candidate.get_alternative_allele());
//    });
//}
//
//HaplotypePhaser::Haplotypes HaplotypePhaser::get_unique_haplotypes_from_tree(const GenomicRegion& a_region)
//{
//    auto haplotypes = the_tree_.get_haplotypes(a_region);
//    unique_least_complex(haplotypes);
//    return haplotypes;
//}
//
////HaplotypePhaser::HaplotypePseudoCounts
////HaplotypePhaser::get_haplotype_prior_counts(const Haplotypes& the_haplotypes,
////                                            CandidateIterator first, CandidateIterator last,
////                                            const GenomicRegion& the_region) const
////{
////    auto haplotype_priors = get_haplotype_prior_probabilities<RealType>(the_haplotypes, first, last);
////    Haplotype reference {the_reference_, the_region};
////    return BayesianGenotypeModel::get_haplotype_prior_pseudo_counts(haplotype_priors, reference,
////                                                                    reference_bias_);
////}
//
//BayesianGenotypeModel::ReadRanges<SampleIdType, HaplotypePhaser::ReadIterator>
//HaplotypePhaser::get_read_iterator_ranges(const GenomicRegion& the_region) const
//{
//    BayesianGenotypeModel::ReadRanges<SampleIdType, ReadIterator> reads {};
//    reads.reserve(the_reads_.size());
//    
//    for (const auto& sample_reads : the_reads_) {
//        reads.emplace(sample_reads.first, overlap_range(sample_reads.second.cbegin(),
//                                                        sample_reads.second.cend(), the_region));
//    }
//    
//    return reads;
//}
//
//void HaplotypePhaser::remove_unlikely_haplotypes_from_tree(const Haplotypes& current_haplotypes,
//                                                           const HaplotypePseudoCounts& prior_counts,
//                                                           const HaplotypePseudoCounts& posterior_counts)
//{
//    for (const auto& haplotype : current_haplotypes) {
//        
////        haplotype.print_explicit_alleles(); cout << endl;
////        cout << prior_counts.at(haplotype) << " " << posterior_counts.at(haplotype) << endl;
//        
//        if (!is_haplotype_likely(BayesianGenotypeModel::haplotype_population_probability(haplotype, prior_counts),
//                                 BayesianGenotypeModel::haplotype_population_probability(haplotype, posterior_counts))) {
//            the_tree_.prune_all(haplotype);
//        } else {
//            // we may as-well clean the tree up here too
//            the_tree_.prune_unique(haplotype);
//        }
//    }
//}
//
//bool HaplotypePhaser::is_haplotype_likely(RealType haplotype_population_prior, RealType haplotype_population_posterior) const
//{
//    if (haplotype_population_posterior < 0.001) return false;
//    if (haplotype_population_posterior > 0.05) return true;
//    
//    auto haplotype_liklihood = haplotype_population_posterior / haplotype_population_prior;
//    
//    if (haplotype_liklihood < 1.0) return false;
//    
//    return true;
//}
//
//void HaplotypePhaser::clear_phased_region(CandidateIterator first, CandidateIterator last)
//{
//    the_candidates_.erase(the_candidates_.cbegin(), last);
//    
//    for (auto& sample_reads : the_reads_) {
//        auto overlapped_reads = overlap_range(sample_reads.second.cbegin(), sample_reads.second.cend(), *first);
//        sample_reads.second.erase(sample_reads.second.cbegin(), overlapped_reads.begin().base());
//    }
//    
//    the_tree_.clear();
//    the_model_.clear_cache();
//}
//
//bool HaplotypePhaser::is_potentially_phasable(const GenomicRegion& current_region, const Variant& previous_candidate) const
//{
//    return has_shared(the_reads_, current_region, previous_candidate);
//}
//
//bool HaplotypePhaser::was_phased(const GenomicRegion& a_region, const HaplotypePseudoCounts& posterior_counts) const
//{
//    return false;
//}

} // namespace Octopus
