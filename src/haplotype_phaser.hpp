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

namespace Octopus
{

class HaplotypePhaser
{
public:
    using SampleGenotypePosteriors = std::unordered_map<Genotype<Haplotype>, double>;
    using GenotypePosteriors       = std::unordered_map<SampleIdType, SampleGenotypePosteriors>;
    
    struct PhaseRegion
    {
        PhaseRegion() = default;
        template <typename Region> PhaseRegion(Region&& region, double probability);
        
        GenomicRegion region;
        double probability;
    };
    using PhaseSet  = std::vector<PhaseRegion>;
    using PhaseSets = std::unordered_map<SampleIdType, PhaseSet>;
    
    HaplotypePhaser() = delete;
    explicit HaplotypePhaser(ReferenceGenome& reference, const std::vector<Variant>& candidates, const ReadMap& reads,
                             unsigned max_haplotypes = 128, unsigned max_indicators = 3);
    ~HaplotypePhaser() = default;
    
    bool done() const noexcept;
    
    std::vector<Haplotype> get_haplotypes();
    void unique(const std::vector<Haplotype>& haplotypes);
    
    PhaseSets phase(const std::vector<Haplotype>& haplotypes, const GenotypePosteriors& genotype_posteriors);
    
private:
    HaplotypeTree tree_;
    
    MappableSet<Variant> buffered_candidates_;
    const ReadMap* reads_;
    
    GenomeWalker walker_;
    
    bool is_phasing_enabled_;
    
    GenomicRegion tree_region_;
    GenomicRegion next_region_;
};

} // namespace Octopus

#endif /* defined(__Octopus__haplotype_phaser__) */
