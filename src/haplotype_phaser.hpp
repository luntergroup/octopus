//
//  haplotype_phaser.hpp
//  Octopus
//
//  Created by Daniel Cooke on 30/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__haplotype_phaser__
#define __Octopus__haplotype_phaser__

#include <cstddef> // std::size_t
#include <string>
#include <vector>
#include <deque>
#include <unordered_map>
#include <boost/range/iterator_range_core.hpp>

#include "common.hpp"
#include "genomic_region.hpp"
#include "map_utils.hpp"
#include "reference_genome.hpp"
#include "aligned_read.hpp"
#include "variant.hpp"
#include "haplotype_tree.hpp"

#include <iostream> // TEST
using std::cout;    // TEST
using std::endl;    // TEST

namespace Octopus
{

//class HaplotypePhaser
//{
//public:
//    using RealType     = Octopus::ProbabilityType;
//    using SampleIdType = Octopus::SampleIdType;
//    template <typename SampleIdType, typename Iterator>
//    using ReadRanges = std::unordered_map<SampleIdType, boost::iterator_range<Iterator>>;
//    
//    struct PhasedRegion
//    {
//        GenomicRegion the_region;
//        std::vector<Haplotype> the_haplotypes;
//        std::vector<Genotype<Haplotype>> the_genotypes;
//        BayesianGenotypeModel::Latents<SampleIdType, RealType> the_latent_posteriors;
//        
//        PhasedRegion() = default;
//        template <typename GenomicRegion_, typename Haplotypes_, typename Genotypes_, typename Latents_>
//        PhasedRegion(GenomicRegion_&& the_region, Haplotypes_&& the_haplotypes,
//                     Genotypes_&& the_genotypes, Latents_&& the_latent_posteriors)
//        :
//        the_region {std::forward<GenomicRegion_>(the_region)},
//        the_haplotypes {std::forward<Haplotypes_>(the_haplotypes)},
//        the_genotypes {std::forward<Genotypes_>(the_genotypes)},
//        the_latent_posteriors {std::forward<Latents_>(the_latent_posteriors)}
//        {}
//    };
//    
//    using PhasedRegions = std::vector<PhasedRegion>;
//    
//    enum class SteamingStatus { Finished, Unfinished };
//    
//    HaplotypePhaser() = delete;
//    HaplotypePhaser(ReferenceGenome& the_reference, VariationalBayesGenotypeModel& the_model,
//                    unsigned ploidy, unsigned max_haplotypes=128, unsigned max_model_update_iterations=3,
//                    RealType reference_bias=1.0);
//    ~HaplotypePhaser() = default;
//    
//    HaplotypePhaser(const HaplotypePhaser&)            = delete;
//    HaplotypePhaser& operator=(const HaplotypePhaser&) = delete;
//    HaplotypePhaser(HaplotypePhaser&&)                 = delete;
//    HaplotypePhaser& operator=(HaplotypePhaser&&)      = delete;
//    
//    template <typename ForwardIterator1, typename ForwardIterator2>
//    void put_data(const ReadRanges<SampleIdType, ForwardIterator1>& the_reads,
//                  ForwardIterator2 first_candidate, ForwardIterator2 last_candidate);
//    
//    PhasedRegions get_phased_regions(SteamingStatus status);
//    
//private:
//    using ReadMap               = std::unordered_map<SampleIdType, std::deque<AlignedRead>>;
//    using Haplotypes            = HaplotypeTree::Haplotypes;
//    using GenotypePosteriors    = BayesianGenotypeModel::GenotypeProbabilities<SampleIdType, RealType>;
//    using HaplotypePseudoCounts = BayesianGenotypeModel::HaplotypePseudoCounts<RealType>;
//    
//    ReadMap the_reads_;
//    std::deque<Variant> the_candidates_;
//    
//    using ReadIterator      = typename ReadMap::mapped_type::const_iterator;
//    using CandidateIterator = decltype(the_candidates_)::const_iterator;
//    
//    ReferenceGenome& the_reference_;
//    unsigned ploidy_;
//    HaplotypeTree the_tree_;
//    VariationalBayesGenotypeModel& the_model_;
//    
//    RealType reference_bias_;
//    
//    unsigned max_haplotypes_;
//    unsigned max_region_density_;
//    unsigned max_model_update_iterations_;
//    
//    PhasedRegion the_last_unphased_region_;
//    PhasedRegions the_phased_regions_;
//    
//    void phase_current_data();
//    void extend_tree(CandidateIterator first, CandidateIterator last);
//    Haplotypes get_unique_haplotypes_from_tree(const GenomicRegion& a_region);
////    HaplotypePseudoCounts get_haplotype_prior_counts(const Haplotypes& the_haplotypes,
////                                                     CandidateIterator first, CandidateIterator last,
////                                                     const GenomicRegion& the_region) const;
//    BayesianGenotypeModel::ReadRanges<SampleIdType, ReadIterator> get_read_iterator_ranges(const GenomicRegion& the_region) const;
//    void remove_unlikely_haplotypes_from_tree(const Haplotypes& current_haplotypes,
//                                              const HaplotypePseudoCounts& prior_counts,
//                                              const HaplotypePseudoCounts& posterior_counts);
//    bool is_haplotype_likely(RealType haplotype_population_prior, RealType haplotype_population_posterior) const;
//    void clear_phased_region(CandidateIterator first, CandidateIterator last);
//    
//    bool is_potentially_phasable(const GenomicRegion& current_region, const Variant& previous_candidate) const;
//    bool was_phased(const GenomicRegion& a_region, const HaplotypePseudoCounts& posterior_counts) const;
//};
//
//template <typename ForwardIterator1, typename ForwardIterator2>
//void HaplotypePhaser::put_data(const ReadRanges<SampleIdType, ForwardIterator1>& the_reads,
//                               ForwardIterator2 first_candidate, ForwardIterator2 last_candidate)
//{
//    for (auto&& map_pair : the_reads) {
//        the_reads_[map_pair.first].insert(the_reads_[map_pair.first].end(), map_pair.second.first, map_pair.second.second);
//    }
//    
//    the_candidates_.insert(the_candidates_.end(), first_candidate, last_candidate);
//    
//    if (the_tree_.empty()) {
//        the_last_unphased_region_.the_region = get_left_overhang(*leftmost_sorted_mappable(the_reads_),
//                                                                 the_candidates_.front().get_region());
//    }
//    
//    phase_current_data();
//}

} // namespace Octopus

#endif /* defined(__Octopus__haplotype_phaser__) */
