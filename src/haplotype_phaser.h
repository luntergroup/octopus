//
//  haplotype_phaser.h
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

#include "genomic_region.h"
#include "map_utils.h"
#include "reference_genome.h"
#include "aligned_read.h"
#include "variant.h"
#include "haplotype_tree.h"
#include "variational_bayes_genotype_model.h"

#include <iostream> // TEST
using std::cout;    // TEST
using std::endl;    // TEST

class HaplotypePhaser
{
public:
    using SampleIdType = std::string;
    template <typename ForwardIterator>
    using ReadRangeMap = std::unordered_map<SampleIdType, std::pair<ForwardIterator, ForwardIterator>>;
    using PhasedGenotypePosteriors = std::vector<VariationalBayesGenotypeModelLatents>;
    
    HaplotypePhaser() = delete;
    HaplotypePhaser(ReferenceGenome& the_reference, VariationalBayesGenotypeModel& the_model,
                    unsigned ploidy);
    ~HaplotypePhaser() = default;
    
    HaplotypePhaser(const HaplotypePhaser&)            = delete;
    HaplotypePhaser& operator=(const HaplotypePhaser&) = delete;
    HaplotypePhaser(HaplotypePhaser&&)                 = delete;
    HaplotypePhaser& operator=(HaplotypePhaser&&)      = delete;
    
    template <typename ForwardIterator1, typename ForwardIterator2>
    void put_data(const ReadRangeMap<ForwardIterator1>& the_reads,
                  ForwardIterator2 first_candidate, ForwardIterator2 last_candidate);
    
    PhasedGenotypePosteriors get_phased_genotypes_posteriors(bool include_partially_phased_haplotypes);
    
private:
    using RealType              = VariationalBayesGenotypeModel::RealType;
    using ReadMap               = std::unordered_map<SampleIdType, std::deque<AlignedRead>>;
    using Haplotypes            = HaplotypeTree::Haplotypes;
    using GenotypePosteriors    = VariationalBayesGenotypeModel::GenotypeResponsabilities;
    using HaplotypePseudoCounts = VariationalBayesGenotypeModel::HaplotypePseudoCounts;
    
    ReadMap the_reads_;
    std::deque<Variant> the_candidates_;
    
    using ReadIterator      = typename ReadMap::mapped_type::const_iterator;
    using CandidateIterator = decltype(the_candidates_)::const_iterator;
    
    ReferenceGenome& the_reference_;
    unsigned ploidy_;
    HaplotypeTree the_tree_;
    VariationalBayesGenotypeModel& the_model_;
    
    std::size_t max_haplotypes_;
    unsigned max_region_density_;
    unsigned max_model_update_iterations_;
    
    PhasedGenotypePosteriors the_phased_genotypes_;
    
    GenomicRegion the_last_unphased_region_;
    VariationalBayesGenotypeModelLatents the_last_unphased_posteriors_;
    
    void phase();
    void extend_haplotypes(CandidateIterator first, CandidateIterator last);
    Haplotypes get_haplotypes(const GenomicRegion& a_region);
    HaplotypePseudoCounts get_haplotype_prior_counts(const Haplotypes& the_haplotypes,
                                                     CandidateIterator first, CandidateIterator last,
                                                     const GenomicRegion& the_region) const;
    SamplesReads<ReadIterator> get_read_ranges(const GenomicRegion& the_region) const;
    void remove_unlikely_haplotypes(const Haplotypes& the_haplotypes,
                                    const HaplotypePseudoCounts& prior_counts,
                                    const HaplotypePseudoCounts& posterior_counts);
    void remove_phased_region(CandidateIterator first, CandidateIterator last);
};

template <typename ForwardIterator1, typename ForwardIterator2>
void HaplotypePhaser::put_data(const ReadRangeMap<ForwardIterator1>& the_reads,
                               ForwardIterator2 first_candidate, ForwardIterator2 last_candidate)
{
    for (auto&& map_pair : the_reads) {
        the_reads_[map_pair.first].insert(the_reads_[map_pair.first].end(), map_pair.second.first, map_pair.second.second);
    }
    
    the_candidates_.insert(the_candidates_.end(), first_candidate, last_candidate);
    
    if (the_tree_.empty()) {
        the_last_unphased_region_ = get_left_overhang(*leftmost_sorted_mappable(the_reads_),
                                                      the_candidates_.front().get_region());
    }
    
    phase();
}

#endif /* defined(__Octopus__haplotype_phaser__) */
