// Copyright (c) 2017 Daniel Cooke
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
#include "basics/genomic_region.hpp"
#include "containers/probability_matrix.hpp"
#include "concepts/mappable.hpp"
#include "concepts/mappable_range.hpp"
#include "utils/mappable_algorithms.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "basics/phred.hpp"

namespace octopus {

class Phaser
{
public:
    using GenotypePosteriorMap       = ProbabilityMatrix<Genotype<Haplotype>>;
    using SampleGenotypePosteriorMap = GenotypePosteriorMap::InnerMap;
    using GenotypeCallMap            = std::unordered_map<SampleName, Genotype<Haplotype>>;
    
    struct PhaseSet;
    
    Phaser() = default;
    
    explicit Phaser(Phred<double> min_phase_score);
    
    Phaser(const Phaser&)            = default;
    Phaser& operator=(const Phaser&) = default;
    Phaser(Phaser&&)                 = default;
    Phaser& operator=(Phaser&&)      = default;
    
    ~Phaser() = default;
    
    boost::optional<PhaseSet> try_phase(const std::vector<Haplotype>& haplotypes,
                                        const GenotypePosteriorMap& genotype_posteriors,
                                        const std::vector<GenomicRegion>& regions) const;
    
    PhaseSet force_phase(const std::vector<Haplotype>& haplotypes,
                         const GenotypePosteriorMap& genotype_posteriors,
                         const std::vector<GenomicRegion>& regions,
                         boost::optional<GenotypeCallMap> genotype_calls = boost::none) const;
    
private:
    Phred<double> min_phase_score_;
    boost::optional<Phred<double>> max_phase_score_ = Phred<double> {100};
};

struct Phaser::PhaseSet
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
