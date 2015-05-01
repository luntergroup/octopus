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

class HaplotypePhaser
{
public:
    using SampleIdType = std::string;
    using PhasedGenotypePosteriors = std::vector<GenotypePosteriors>;
    
    HaplotypePhaser() = delete;
    HaplotypePhaser(ReferenceGenome& the_reference, VariationalBayesGenotypeModel& the_model,
                    unsigned ploidy);
    ~HaplotypePhaser() = default;
    
    HaplotypePhaser(const HaplotypePhaser&)            = delete;
    HaplotypePhaser& operator=(const HaplotypePhaser&) = delete;
    HaplotypePhaser(HaplotypePhaser&&)                 = delete;
    HaplotypePhaser& operator=(HaplotypePhaser&&)      = delete;
    
    template <typename ForwardIterator1, typename ForwardIterator2>
    void put_data(SampleIdType sample, ForwardIterator1 first_read, ForwardIterator1 last_read,
                  ForwardIterator2 first_candidate, ForwardIterator2 last_candidate);
    
    PhasedGenotypePosteriors get_phased_genotypes_posteriors(bool include_partially_phased_haplotypes);
    
private:
    using RealType          = VariationalBayesGenotypeModel::RealType;
    using ReadMap           = std::unordered_map<SampleIdType, std::deque<AlignedRead>>;
    
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
    unsigned max_indicators_;
    unsigned max_model_update_iterations_;
    
    PhasedGenotypePosteriors the_phased_genotypes_;
    
    GenomicRegion the_last_unphased_region_;
    
    void phase();
    void remove_phased_region(CandidateIterator first, CandidateIterator last);
};

template <typename ForwardIterator1, typename ForwardIterator2>
void HaplotypePhaser::put_data(SampleIdType sample, ForwardIterator1 first_read, ForwardIterator1 last_read,
                               ForwardIterator2 first_candidate, ForwardIterator2 last_candidate)
{
    the_reads_[sample].insert(the_reads_[sample].end(), first_read, last_read);
    the_candidates_.insert(the_candidates_.end(), first_candidate, last_candidate);
    
//    if (the_tree_.empty()) {
//        the_last_unphased_region_ = leftmost_sorted_mappable(the_candidates_)->get_region();
//    }
    
    phase();
}

#endif /* defined(__Octopus__haplotype_phaser__) */
