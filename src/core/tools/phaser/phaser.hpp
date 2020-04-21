// Copyright (c) 2015-2020 Daniel Cooke
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
#include "core/types/indexed_haplotype.hpp"
#include "core/types/genotype.hpp"
#include "utils/mappable_algorithms.hpp"

namespace octopus {

class Phaser
{
public:
    using GenotypePosteriorMap       = ProbabilityMatrix<Genotype<IndexedHaplotype<>>>;
    using SampleGenotypePosteriorMap = GenotypePosteriorMap::InnerMap;
    using GenotypeCallMap            = std::unordered_map<SampleName, Genotype<IndexedHaplotype<>>>;
    
    enum class GenotypeMatchType { exact, unique };
    
    struct Config
    {
        GenotypeMatchType genotype_match = GenotypeMatchType::exact;
        Phred<double> min_phase_quality = Phred<double> {10};
        boost::optional<Phred<double>> max_phase_quality = Phred<double> {100};
    };
    
    struct PhaseSet
    {
        std::vector<std::size_t> site_indices;
        Phred<double> quality;
    };
    using PhaseSetVector = std::vector<PhaseSet>;
    using PhaseSetMap = std::unordered_map<SampleName, PhaseSetVector>;
    
    Phaser() = default;
    
    Phaser(Config config);
    
    Phaser(const Phaser&)            = default;
    Phaser& operator=(const Phaser&) = default;
    Phaser(Phaser&&)                 = default;
    Phaser& operator=(Phaser&&)      = default;
    
    ~Phaser() = default;
    
    PhaseSetMap
    phase(const MappableBlock<Haplotype>& haplotypes,
          const GenotypePosteriorMap& genotype_posteriors,
          const std::vector<GenomicRegion>& variation_regions,
          boost::optional<GenotypeCallMap> genotype_calls = boost::none) const;
    
private:
    using CompressedGenotype = Genotype<IndexedHaplotype<>>;
    
    Config config_;
    
    PhaseSetVector
    phase_sample(const std::vector<GenomicRegion>& partitions,
                 const std::vector<CompressedGenotype>& genotypes,
                 const SampleGenotypePosteriorMap& genotype_posteriors) const;
};

namespace debug {

template <typename S>
void print_phase_sets(S&& stream, const Phaser::PhaseSetMap& phasings, const std::vector<GenomicRegion>& variation_regions)
{
    if (!variation_regions.empty()) {
        stream << "Phase sets in " << encompassing_region(variation_regions) << '\n';
        for (const auto& p : phasings) {
            if (!p.second.empty()) {
                stream << "Phase sets for sample " << p.first << '\n';
                for (const auto& phase_set : p.second) {
                    stream << "\t*";
                    for (auto site_idx : phase_set.site_indices) {
                        stream << variation_regions[site_idx] << " ";
                    }
                    stream << phase_set.quality << '\n';
                }
            }
        }
    }
}

void print_phase_sets(const Phaser::PhaseSetMap& phasings, const std::vector<GenomicRegion>& variation_regions);

} // namespace debug

} // namespace octopus

#endif
