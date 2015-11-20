//
//  haplotype_phaser.hpp
//  Octopus
//
//  Created by Daniel Cooke on 30/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__haplotype_phaser__
#define __Octopus__haplotype_phaser__

#include <vector>
#include <unordered_map>

#include "common.hpp"
#include "genomic_region.hpp"
#include "variant.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "reference_genome.hpp"
#include "haplotype_tree.hpp"
#include "mappable_set.hpp"
#include "genome_walker.hpp"

//namespace Octopus
//{

using namespace Octopus; // for now

class HaplotypePhaser
{
public:
    using UnphasedSampleGenotypePosteriors = std::unordered_map<Genotype<Haplotype>, double>;
    using UnphasedGenotypePosteriors       = std::unordered_map<SampleIdType, UnphasedSampleGenotypePosteriors>;
    using PhasedSampleGenotypePosteriors   = std::unordered_map<Genotype<Haplotype>, double>;
    using PhasedGenotypePosteriors         = std::unordered_map<SampleIdType, PhasedSampleGenotypePosteriors>;
    
    HaplotypePhaser() = delete;
    HaplotypePhaser(ReferenceGenome& reference);
    HaplotypePhaser(ReferenceGenome& reference, unsigned max_haplotypes);
    HaplotypePhaser(ReferenceGenome& reference, unsigned max_haplotypes, unsigned max_indicators);
    
    void setup(const std::vector<Variant>& candidates, const ReadMap& reads);
    bool expended_candidates() const noexcept;
    std::vector<Haplotype> get_haplotypes() const;
    PhasedGenotypePosteriors phase(const std::vector<Haplotype>& haplotypes,
                                   const UnphasedGenotypePosteriors& genotype_posteriors,
                                   const ReadMap& reads);
    
private:
    HaplotypeTree tree_;
    MappableSet<Variant> buffered_candidates_;
    GenomicRegion tree_region_; // until HaplotypeTree supports get_haplotypes()
    
    unsigned max_haplotypes_ = 10'000;
    unsigned max_indicators_ = 5;
    
    GenomeWalker walker_;
    
    void extend_tree(const ReadMap& reads);
};

//} // namespace Octopus

#endif /* defined(__Octopus__haplotype_phaser__) */
