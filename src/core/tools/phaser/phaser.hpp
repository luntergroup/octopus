// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef phaser_hpp
#define phaser_hpp

#include <vector>
#include <unordered_map>
#include <utility>
#include <functional>
#include <cassert>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "concepts/mappable.hpp"
#include "concepts/mappable_range.hpp"
#include "basics/genomic_region.hpp"
#include "basics/phred.hpp"
#include "containers/probability_matrix.hpp"
#include "containers/mappable_block.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "utils/mappable_algorithms.hpp"

namespace octopus {

class Phaser
{
public:
    using GenotypePosteriorMap       = ProbabilityMatrix<Genotype<Haplotype>>;
    using SampleGenotypePosteriorMap = GenotypePosteriorMap::InnerMap;
    using GenotypeCallMap            = std::unordered_map<SampleName, Genotype<Haplotype>>;
    
    enum class GenotypeMatchType { exact, unique };
    
    struct Config
    {
        GenotypeMatchType genotype_match = GenotypeMatchType::exact;
        Phred<double> min_phase_score = Phred<double> {10};
        boost::optional<Phred<double>> max_phase_score = Phred<double> {100};
    };
    
    struct PhaseSet
    {
        struct PhaseRegion : public Mappable<PhaseRegion>
        {
            PhaseRegion() = default;
            
            template <typename Region> PhaseRegion(Region&& region, Phred<double> score);
            
            GenomicRegion region;
            Phred<double> score;
            
            const GenomicRegion& mapped_region() const noexcept { return region; }
        };
        
        using SamplePhaseRegions = std::vector<PhaseRegion>;
        using PhaseRegions       = std::unordered_map<SampleName, SamplePhaseRegions>;
        
        PhaseSet() = delete;
        
        template <typename R> PhaseSet(R&& region);
        template <typename R, typename T> PhaseSet(R&& region, T&& phase_regions);
        
        GenomicRegion region;
        PhaseRegions phase_regions;
    };
    
    Phaser() = default;
    
    Phaser(Config config);
    
    Phaser(const Phaser&)            = default;
    Phaser& operator=(const Phaser&) = default;
    Phaser(Phaser&&)                 = default;
    Phaser& operator=(Phaser&&)      = default;
    
    ~Phaser() = default;
    
    PhaseSet
    phase(const MappableBlock<Haplotype>& haplotypes,
          const GenotypePosteriorMap& genotype_posteriors,
          const std::vector<GenomicRegion>& variation_regions,
          boost::optional<GenotypeCallMap> genotype_calls = boost::none) const;
    
private:
    using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;
    
    Config config_;
    
    PhaseSet::SamplePhaseRegions
    phase_sample(const GenomicRegion& region,
                 const std::vector<GenomicRegion>& partitions,
                 const std::vector<GenotypeReference>& genotypes,
                 const SampleGenotypePosteriorMap& genotype_posteriors) const;
};

template <typename Region>
Phaser::PhaseSet::PhaseRegion::PhaseRegion(Region&& region, const Phred<double> score)
: region {std::forward<Region>(region)}
, score {score}
{}

template <typename R>
Phaser::PhaseSet::PhaseSet(R&& region)
: region {std::forward<R>(region)}
, phase_regions {}
{}

template <typename R, typename T>
Phaser::PhaseSet::PhaseSet(R&& region, T&& phase_regions)
: region {std::forward<R>(region)}
, phase_regions {std::forward<T>(phase_regions)}
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

bool is_split_phasing(const Phaser::PhaseSet& phase);

namespace debug {

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

} // namespace debug

} // namespace octopus

#endif
