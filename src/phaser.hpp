//
//  phaser.hpp
//  Octopus
//
//  Created by Daniel Cooke on 20/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef phaser_hpp
#define phaser_hpp

#include <vector>
#include <unordered_map>
#include <utility>

#include <boost/optional.hpp>

#include "common.hpp"
#include "genomic_region.hpp"
#include "probability_matrix.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"

namespace Octopus
{
    class Phaser
    {
    public:
        using GenotypePosteriorMap       = ProbabilityMatrix<Genotype<Haplotype>>;
        using SampleGenotypePosteriorMap = GenotypePosteriorMap::InnerMap;
        
        Phaser() = default;
        
        explicit Phaser(double min_phase_score);
        
        ~Phaser() = default;
        
        Phaser(const Phaser&)            = default;
        Phaser& operator=(const Phaser&) = default;
        Phaser(Phaser&&)                 = default;
        Phaser& operator=(Phaser&&)      = default;
        
        struct PhaseSet
        {
            struct PhaseRegion
            {
                PhaseRegion() = default;
                template <typename Region> PhaseRegion(Region&& region, double score);
                ~PhaseRegion() = default;
                
                PhaseRegion(const PhaseRegion&)            = default;
                PhaseRegion& operator=(const PhaseRegion&) = default;
                PhaseRegion(PhaseRegion&&)                 = default;
                PhaseRegion& operator=(PhaseRegion&&)      = default;
                
                GenomicRegion region;
                double score;
            };
            
            using SamplePhaseRegions = std::vector<PhaseRegion>;
            using PhaseRegions       = std::unordered_map<SampleIdType, SamplePhaseRegions>;
            
            PhaseSet() = delete;
            template <typename R> PhaseSet(R&& region);
            template <typename R, typename T> PhaseSet(R&& region, T&& phase_regions);
            ~PhaseSet() = default;
            
            PhaseSet(const PhaseSet&)            = default;
            PhaseSet& operator=(const PhaseSet&) = default;
            PhaseSet(PhaseSet&&)                 = default;
            PhaseSet& operator=(PhaseSet&&)      = default;
            
            GenomicRegion region;
            PhaseRegions phase_regions;
        };
        
        boost::optional<PhaseSet> try_phase(const std::vector<Haplotype>& haplotypes,
                                            const GenotypePosteriorMap& genotype_posteriors);
        
        PhaseSet force_phase(const std::vector<Haplotype>& haplotypes,
                             const GenotypePosteriorMap& genotype_posteriors);
        
    private:
        double min_phase_score_;
    };
    
    template <typename Region>
    Phaser::PhaseSet::PhaseRegion::PhaseRegion(Region&& region, const double score)
    : region {std::forward<Region>(region)}, score {score} {}
    
    template <typename R>
    Phaser::PhaseSet::PhaseSet(R&& region)
    : region {std::forward<R>(region)}, phase_regions {} {}
    
    template <typename R, typename T>
    Phaser::PhaseSet::PhaseSet(R&& region, T&& phase_regions)
    : region {std::forward<R>(region)}, phase_regions {std::forward<T>(phase_regions)} {}
} // namespace Octopus

#endif /* phaser_hpp */
