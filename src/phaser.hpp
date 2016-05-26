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
#include <functional>
#include <cassert>

#include <boost/optional.hpp>

#include "common.hpp"
#include "genomic_region.hpp"
#include "probability_matrix.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "variant.hpp"
#include "mappable.hpp"
#include "mappable_ranges.hpp"
#include "mappable_algorithms.hpp"

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
            struct PhaseRegion : public Mappable<PhaseRegion>
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
                
                const GenomicRegion& mapped_region() const noexcept { return region; }
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
                                            const GenotypePosteriorMap& genotype_posteriors,
                                            const std::vector<Variant>& candidates);
        
        PhaseSet force_phase(const std::vector<Haplotype>& haplotypes,
                             const GenotypePosteriorMap& genotype_posteriors,
                             const std::vector<Variant>& candidates);
        
    private:
        double min_phase_score_;
    };
    
    template <typename Region>
    Phaser::PhaseSet::PhaseRegion::PhaseRegion(Region&& region, const double score)
    :
    region {std::forward<Region>(region)},
    score {score}
    {}
    
    template <typename R>
    Phaser::PhaseSet::PhaseSet(R&& region)
    :
    region {std::forward<R>(region)},
    phase_regions {}
    {}
    
    template <typename R, typename T>
    Phaser::PhaseSet::PhaseSet(R&& region, T&& phase_regions)
    :
    region {std::forward<R>(region)},
    phase_regions {std::forward<T>(phase_regions)}
    {}
    
    // non-member methods
    
    template <typename T>
    boost::optional<std::reference_wrapper<const Phaser::PhaseSet::PhaseRegion>>
    find_phase_region(const Phaser::PhaseSet::SamplePhaseRegions& phasings, const T& mappable)
    {
        const auto overlapped = overlap_range(phasings, mappable, BidirectionallySortedTag {});
        
        if (!overlapped.empty()) {
            assert(size(overlapped, BidirectionallySortedTag {}) == 1);
            return std::cref(overlapped.front());
        }
        
        return boost::none;
    }
    
    namespace debug
    {
        template <typename S>
        void print_phase_sets(S&& stream, const Phaser::PhaseSet& phasings)
        {
            stream << "Phased region is " << phasings.region << '\n';
            
            for (const auto& p : phasings.phase_regions) {
                if (p.second.size() > 1) {
                    stream << "Phase regions for sample " << p.first << '\n';
                    for (const auto& r : p.second) {
                        stream << "\t* " << r.region << " " << r.score << '\n';
                    }
                }
            }
        }
        
        void print_phase_sets(const Phaser::PhaseSet& phasings);
    }
} // namespace Octopus

#endif /* phaser_hpp */
