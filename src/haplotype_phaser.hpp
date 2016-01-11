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
#include <utility>

#include <boost/optional.hpp>

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
    
    struct PhaseSet
    {
        struct PhaseRegion
        {
            PhaseRegion() = default;
            template <typename Region> PhaseRegion(Region&& region, double score);
            
            GenomicRegion region;
            double score;
        };
        
        using SamplePhaseRegions = std::vector<PhaseRegion>;
        using PhaseRegions       = std::unordered_map<SampleIdType, SamplePhaseRegions>;
        
        PhaseSet() = default;
        template <typename R> PhaseSet(R&& region);
        template <typename R, typename T> PhaseSet(R&& region, T&& phase_regions);
        
        GenomicRegion region;
        PhaseRegions phase_regions;
    };
    
    HaplotypePhaser() = delete;
    
    explicit HaplotypePhaser(const ReferenceGenome& reference,
                             const std::vector<Variant>& candidates,
                             const ReadMap& reads,
                             unsigned max_haplotypes = 128,
                             unsigned max_indicators = 3);
    
    ~HaplotypePhaser() = default;
    
    HaplotypePhaser(const HaplotypePhaser&)            = default;
    HaplotypePhaser& operator=(const HaplotypePhaser&) = default;
    HaplotypePhaser(HaplotypePhaser&&)                 = default;
    HaplotypePhaser& operator=(HaplotypePhaser&&)      = default;
    
    bool done() const noexcept;
    
    std::vector<Haplotype> get_haplotypes();
    std::vector<Haplotype> get_haplotypes(const GenotypePosteriors& genotype_posteriors);
    
    void set_haplotypes(const std::vector<Haplotype>& haplotypes);
    
    boost::optional<PhaseSet> phase(const std::vector<Haplotype>& haplotypes,
                                    const GenotypePosteriors& genotype_posteriors);
    
private:
    HaplotypeTree tree_;
    
    MappableSet<Variant> buffered_candidates_;
    const ReadMap* reads_;
    
    GenomeWalker walker_;
    
    bool is_phasing_enabled_;
    
    GenomicRegion current_region_, next_region_;
    
    void remove_low_posterior_haplotypes(const std::vector<Haplotype>& haplotypes,
                                         const GenotypePosteriors& genotype_posteriors);
    
    PhaseSet find_optimal_phase_set(const GenomicRegion& region,
                                    MappableSet<Variant> variants,
                                    const GenotypePosteriors& genotype_posteriors);
};

template <typename Region>
HaplotypePhaser::PhaseSet::PhaseRegion::PhaseRegion(Region&& region, double score)
: region {std::forward<Region>(region)}, score {score} {}

template <typename R>
HaplotypePhaser::PhaseSet::PhaseSet(R&& region)
: region {std::forward<R>(region)}, phase_regions {} {}

template <typename R, typename T>
HaplotypePhaser::PhaseSet::PhaseSet(R&& region, T&& phase_regions)
: region {std::forward<R>(region)}, phase_regions {std::forward<T>(phase_regions)} {}

} // namespace Octopus

#endif /* defined(__Octopus__haplotype_phaser__) */
