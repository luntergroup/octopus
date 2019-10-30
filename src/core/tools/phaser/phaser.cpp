// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "phaser.hpp"

#include <deque>
#include <algorithm>
#include <numeric>
#include <cstddef>
#include <cmath>
#include <utility>
#include <iostream>

#include <boost/functional/hash.hpp>

#include "utils/mappable_algorithms.hpp"
#include "utils/maths.hpp"

#include "timers.hpp"

namespace octopus {

Phaser::Phaser(Config config) : config_ {config} {}

namespace {

using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;

std::vector<GenotypeReference> extract_genotypes(const Phaser::GenotypePosteriorMap& genotype_posteriors)
{
    return extract_key_refs(genotype_posteriors);
}

auto minmax_ploidy(const std::vector<GenotypeReference>& genotypes) noexcept
{
    assert(!genotypes.empty());
    const static auto ploidy_less = [] (const auto& lhs, const auto& rhs) { return lhs.get().ploidy() < rhs.get().ploidy(); };
    auto p = std::minmax_element(std::cbegin(genotypes), std::cend(genotypes), ploidy_less);
    return std::make_pair(p.first->get().ploidy(), p.second->get().ploidy());
}

double maximum_entropy(const std::size_t num_elements)
{
    return std::log2(num_elements);
}

auto min_phase_score(double p)
{
    if (maths::almost_one(p)) {
        static const Phred<double> max_possible_score {Phred<double>::Probability {0.0}};
        return max_possible_score;
    } else {
        p = std::max(p, 0.5);
        const auto e = -(p * std::log2(p) + (1.0 - p) * std::log2(1.0 - p)) / maximum_entropy(2);
        return Phred<double> {Phred<double>::Probability {std::max(e, 0.0)}};
    }
}

auto min_phase_score(const Genotype<Haplotype>& called_genotype, const Phaser::SampleGenotypePosteriorMap& genotype_posteriors)
{
    return min_phase_score(genotype_posteriors[called_genotype]);
}

Phaser::GenotypePosteriorMap marginalise_collapsed_genotypes(const Phaser::GenotypePosteriorMap& genotype_posteriors)
{
    auto genotypes = extract_key_refs(genotype_posteriors);
    const auto num_genotypes = genotypes.size();
    std::vector<Genotype<Haplotype>> collapsed_genotypes {};
    collapsed_genotypes.reserve(num_genotypes);
    std::unordered_map<Genotype<Haplotype>, std::size_t> collapsed_genotype_indices {};
    collapsed_genotype_indices.reserve(num_genotypes);
    std::vector<std::size_t> collapsed_genotype_index_table {};
    collapsed_genotype_index_table.reserve(num_genotypes);
    for (const auto& genotype : genotypes) {
        auto collapsed = genotype.get().collapse();
        const auto itr = collapsed_genotype_indices.find(collapsed);
        if (itr == std::cend(collapsed_genotype_indices)) {
            const auto index = collapsed_genotypes.size();
            collapsed_genotype_indices.emplace(collapsed, index);
            collapsed_genotypes.push_back(std::move(collapsed));
            collapsed_genotype_index_table.push_back(index);
        } else {
            collapsed_genotype_index_table.push_back(itr->second);
        }
    }
    genotypes.clear();
    genotypes.shrink_to_fit();
    collapsed_genotype_indices.clear();
    const auto num_collapsed_genotypes = collapsed_genotypes.size();
    Phaser::GenotypePosteriorMap result {std::make_move_iterator(std::begin(collapsed_genotypes)),
                                         std::make_move_iterator(std::end(collapsed_genotypes))};
    collapsed_genotypes.clear();
    collapsed_genotypes.shrink_to_fit();
    for (const auto& sample_posteriors : genotype_posteriors) {
        std::vector<double> posteriors(num_collapsed_genotypes);
        std::size_t genotype_idx {0};
        std::for_each(std::cbegin(sample_posteriors.second), std::cend(sample_posteriors.second), [&] (const auto& p) {
            posteriors[collapsed_genotype_index_table[genotype_idx]] += p.second;
            ++genotype_idx;
        });
        insert_sample(sample_posteriors.first, std::move(posteriors), result);
    }
    return result;
}

} // namespace

Phaser::PhaseSet
Phaser::phase(const MappableBlock<Haplotype>& haplotypes,
              const GenotypePosteriorMap& genotype_posteriors,
              const std::vector<GenomicRegion>& variation_regions,
              boost::optional<GenotypeCallMap> genotype_calls) const
{
    assert(!haplotypes.empty());
    assert(!genotype_posteriors.empty1() && !genotype_posteriors.empty2());
    assert(std::is_sorted(std::cbegin(variation_regions), std::cend(variation_regions)));
    const auto& haplotype_region = mapped_region(haplotypes);
    const auto partitions = extract_covered_regions(variation_regions);
    auto genotypes = extract_genotypes(genotype_posteriors);
    unsigned min_genotype_ploidy, max_genotype_ploidy; std::tie(min_genotype_ploidy, max_genotype_ploidy) = minmax_ploidy(genotypes);
    PhaseSet result {haplotype_region};
    result.phase_regions.reserve(genotype_posteriors.size1());
    if (max_genotype_ploidy == 1 || partitions.size() == 1) {
        for (const auto& p : genotype_posteriors) {
            if (config_.max_phase_score) {
                result.phase_regions[p.first].emplace_back(haplotype_region, *config_.max_phase_score);
            } else {
                static const Phred<double> max_possible_score {Phred<double>::Probability {0.0}};
                result.phase_regions[p.first].emplace_back(haplotype_region, max_possible_score);
            }
        }
    } else {
        boost::optional<GenotypePosteriorMap> collapsed_genotype_posteriors {};
        for (const auto& p : genotype_posteriors) {
            const SampleName& sample {p.first};
            if (genotype_calls && config_.max_phase_score && min_phase_score(genotype_calls->at(sample), p.second) >= *config_.max_phase_score) {
                result.phase_regions[sample].emplace_back(haplotype_region, *config_.max_phase_score);
            } else {
                if (!collapsed_genotype_posteriors && (max_genotype_ploidy > 2 || min_genotype_ploidy != max_genotype_ploidy)) {
                    collapsed_genotype_posteriors = marginalise_collapsed_genotypes(genotype_posteriors);
                    genotypes = extract_genotypes(*collapsed_genotype_posteriors);
                    std::tie(min_genotype_ploidy, max_genotype_ploidy) = minmax_ploidy(genotypes);
                }
                PhaseSet::SamplePhaseRegions phases;
                if (collapsed_genotype_posteriors) {
                    phases = phase_sample(haplotype_region, partitions, genotypes, (*collapsed_genotype_posteriors)[sample]);
                } else {
                    phases = phase_sample(haplotype_region, partitions, genotypes, p.second);
                }
                if (config_.max_phase_score) {
                    for (auto& phase : phases) phase.score = std::min(phase.score, *config_.max_phase_score);
                }
                result.phase_regions.emplace(sample, std::move(phases));
            }
        }
    }
    return result;
}

namespace {

using PhaseComplementSet   = std::deque<GenotypeReference>;
using PhaseComplementSets  = std::vector<PhaseComplementSet>;
using PartitionIterator    = std::vector<GenomicRegion>::const_iterator;
using GenotypeChunk        = Genotype<Haplotype>;
using GenotypeChunkVector  = std::vector<GenotypeChunk>;

struct GenotypeChunkPartitions
{
    GenotypeReference genotype;
    GenotypeChunkVector chunks;
};

using GenotypeChunkTable = std::vector<GenotypeChunkPartitions>;

template <typename Container>
GenotypeChunkTable
make_chunk_table(const Container& genotypes, PartitionIterator first_parition, PartitionIterator last_partition)
{
    const auto num_partitions = std::distance(first_parition, last_partition);
    GenotypeChunkTable result {};
    result.reserve(genotypes.size());
    for (const auto& genotype : genotypes) {
        result.push_back({std::cref(genotype), GenotypeChunkVector {}});
        result.back().chunks.reserve(num_partitions);
    }
    std::for_each(first_parition, last_partition, [&] (const auto& region) {
        auto chunks = copy_each<Haplotype>(genotypes, region);
        for (std::size_t i {0}; i < genotypes.size(); ++i) {
            result[i].chunks.push_back(std::move(chunks[i]));
        }
    });
    return result;
}

auto max_ploidy(const GenotypeChunkTable& table) noexcept
{
    assert(!table.empty());
    const static auto ploidy_less = [] (const auto& lhs, const auto& rhs) { return lhs.genotype.get().ploidy() < rhs.genotype.get().ploidy(); };
    return std::max_element(std::cbegin(table), std::cend(table), ploidy_less)->genotype.get().ploidy();
}

struct GenotypeChunkVectorExactHash
{
    auto operator()(const GenotypeChunkVector& chunks) const noexcept
    {
        std::size_t result {};
        for (const auto& genotype : chunks) {
            using boost::hash_combine;
            hash_combine(result, std::hash<Genotype<Haplotype>>{}(genotype));
        }
        return result;
    }
};
struct GenotypeUniqueHash
{
    std::size_t operator()(const Genotype<Haplotype>& genotype) const
    {
        std::size_t result {0};
        for (const Haplotype& haplotype : genotype.copy_unique_ref()) {
            using boost::hash_combine;
            hash_combine(result, std::hash<Haplotype>{}(haplotype));
        }
        return result;
    }
};
struct GenotypeChunkVectorUniqueHash
{
    auto operator()(const GenotypeChunkVector& chunks) const noexcept
    {
        std::size_t result {};
        for (const auto& genotype : chunks) {
            using boost::hash_combine;
            hash_combine(result, GenotypeUniqueHash{}(genotype));
        }
        return result;
    }
};
struct GenotypeChunkVectorExactEqual
{
    bool operator()(const GenotypeChunkVector& lhs, const GenotypeChunkVector& rhs) const
    {
        return lhs == rhs;
    }
};
bool have_same_haplotypes(const Genotype<Haplotype>& lhs, const Genotype<Haplotype>& rhs) noexcept
{
    if (lhs.ploidy() < 3 && rhs.ploidy() < 3) {
        return lhs == rhs;
    } else {
        return lhs.copy_unique_ref() == rhs.copy_unique_ref();
    }
}

struct GenotypeChunkUniqueEqual
{
    bool operator()(const GenotypeChunk& lhs, const GenotypeChunk& rhs) const
    {
        return have_same_haplotypes(lhs, rhs);
    }
};
struct GenotypeChunkVectorUniqueEqual
{
    bool operator()(const GenotypeChunkVector& lhs, const GenotypeChunkVector& rhs) const
    {
        return lhs.size() == rhs.size() && std::equal(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs), GenotypeChunkUniqueEqual {});
    }
};

template <typename Hash, typename Equal>
using PhaseSetMap = std::unordered_map<GenotypeChunkVector, PhaseComplementSet, Hash, Equal>;

template <typename Map>
PhaseComplementSets
make_phase_sets_helper(const GenotypeChunkTable& table, Map& phase_sets)
{
    phase_sets.reserve(table.size());
    for (const auto& row : table) {
        phase_sets[row.chunks].push_back(row.genotype);
    }
    PhaseComplementSets result {};
    result.reserve(phase_sets.size());
    for (auto&& p : phase_sets) {
        result.push_back(std::move(p.second));
    }
    return result;
}

PhaseComplementSets
make_phase_sets(const GenotypeChunkTable& table, const Phaser::GenotypeMatchType match_type)
{
    if (match_type == Phaser::GenotypeMatchType::exact || max_ploidy(table) < 3) {
        // if haploid or diploid genotypes then always use exact match as exact and unique match are
        // identical and exact match is faster
        PhaseSetMap<GenotypeChunkVectorExactHash, GenotypeChunkVectorExactEqual> phase_sets {};
        return make_phase_sets_helper(table, phase_sets);
    } else {
        PhaseSetMap<GenotypeChunkVectorUniqueHash, GenotypeChunkVectorUniqueEqual> phase_sets {};
        return make_phase_sets_helper(table, phase_sets);
    }
}

template <typename Container>
PhaseComplementSets
generate_phase_complement_sets(const Container& genotypes,
                               PartitionIterator first_parition, PartitionIterator last_partition,
                               const Phaser::GenotypeMatchType match_type)
{
    const auto table = make_chunk_table(genotypes, first_parition, last_partition);
    return make_phase_sets(table, match_type);
}

template <typename Map>
double marginalise(const PhaseComplementSet& phase_set, const Map& genotype_posteriors)
{
    return std::accumulate(std::cbegin(phase_set), std::cend(phase_set), 0.0,
                           [&] (const auto curr, const auto& genotype) {
                               return curr + genotype_posteriors.at(genotype);
                           });
}

template <typename Map>
double calculate_entropy(const PhaseComplementSet& phase_set, const Map& genotype_posteriors)
{
    const auto norm = marginalise(phase_set, genotype_posteriors);
    if (norm <= 0.0) {
        // if norm ~= 0 then every element in the phase must must have probability ~= 0, so it
        // just looks like a uniform distirbution
        return maximum_entropy(phase_set.size());
    }
    return std::max(0.0, -std::accumulate(std::cbegin(phase_set), std::cend(phase_set), 0.0,
                                          [&genotype_posteriors, norm] (const auto curr, const auto& genotype) {
                                              const auto p = genotype_posteriors.at(genotype) / norm;
                                              return curr + p * std::log2(std::max(p, std::numeric_limits<double>::min()));
                                          }));
}

template <typename Map>
double calculate_relative_entropy(const PhaseComplementSet& phase_set, const Map& genotype_posteriors)
{
    if (phase_set.size() < 2) return 1.0;
    return 1.0 - std::min(calculate_entropy(phase_set, genotype_posteriors) / maximum_entropy(2), 1.0);
}

template <typename Map>
auto calculate_phase_score(const PhaseComplementSet& phase_set, const Map& genotype_posteriors)
{
    return marginalise(phase_set, genotype_posteriors) * calculate_relative_entropy(phase_set, genotype_posteriors);
}

template <typename Map>
Phred<double> calculate_phase_score(const PhaseComplementSets& phase_sets, const Map& genotype_posteriors)
{
    return Phred<double> { Phred<double>::Probability {
    std::max(0.0, 1.0 - std::accumulate(std::cbegin(phase_sets), std::cend(phase_sets), 0.0,
                                        [&] (const auto curr, const auto& phase_set) {
                                            return curr + calculate_phase_score(phase_set, genotype_posteriors);
                                        }))
    }};
}

using GenotypeChunkPosteriorMap = std::unordered_map<Genotype<Haplotype>, double>;

template <typename Container>
auto copy_and_marginalise(const Container& genotypes,
                          const Phaser::SampleGenotypePosteriorMap& genotype_posteriors,
                          const GenomicRegion& region)
{
    auto chunks = copy_each<Haplotype>(genotypes, region);
    GenotypeChunkPosteriorMap chunk_posteriors {genotypes.size()};
    for (std::size_t g {0}; g < genotypes.size(); ++g) {
        chunk_posteriors[chunks[g]] += genotype_posteriors[genotypes[g]];
    }
    std::sort(std::begin(chunks), std::end(chunks), GenotypeLess {});
    chunks.erase(std::unique(std::begin(chunks), std::end(chunks)), std::end(chunks));
    return std::make_pair(std::move(chunks), std::move(chunk_posteriors));
}

template <typename Map>
void print(const PhaseComplementSets& phase_sets, const Map& genotype_posteriors)
{
    for (std::size_t i {0}; i < phase_sets.size(); ++i) {
        const auto& phase_set = phase_sets[i];
        std::cout << "Phase set " << i << ":" << std::endl;;
        for (const Genotype<Haplotype>& genotype : phase_set) {
            octopus::debug::print_variant_alleles(genotype); std::cout << ' ' << genotype_posteriors.at(genotype) << std::endl;
        }
        std::cout << "entropy = " << calculate_entropy(phase_set, genotype_posteriors) << std::endl;
        std::cout << "probability mass = " << marginalise(phase_set, genotype_posteriors) << std::endl;
    }
}

} // namespace

Phaser::PhaseSet::SamplePhaseRegions
Phaser::phase_sample(const GenomicRegion& region,
                     const std::vector<GenomicRegion>& partitions,
                     const std::vector<GenotypeReference>& genotypes,
                     const SampleGenotypePosteriorMap& genotype_posteriors) const
{
    auto first_partition = std::cbegin(partitions);
    auto last_partition  = std::cend(partitions);
    Phaser::PhaseSet::SamplePhaseRegions result {};
    std::vector<Genotype<Haplotype>> chunks;
    GenotypeChunkPosteriorMap chunk_posteriors;
    while (first_partition != std::cend(partitions)) {
        auto curr_region = encompassing_region(first_partition, last_partition);
        std::tie(chunks, chunk_posteriors) = copy_and_marginalise(genotypes, genotype_posteriors, curr_region);
        auto phase_set = generate_phase_complement_sets(chunks, first_partition, last_partition, config_.genotype_match);
        auto phase_score = calculate_phase_score(phase_set, chunk_posteriors);
        if (phase_score >= config_.min_phase_score || std::distance(first_partition, last_partition) == 1) {
            result.emplace_back(encompassing_region(first_partition, last_partition), phase_score);
            first_partition = last_partition;
            last_partition  = std::cend(partitions);
        } else {
            --last_partition;
        }
    }
    return result;
}

// non-member methods

bool is_split_phasing(const Phaser::PhaseSet& phase)
{
    return std::any_of(std::cbegin(phase.phase_regions), std::cend(phase.phase_regions),
                       [] (const auto& p) { return p.second.size() > 1; });
}
    
namespace debug {

void print_phase_sets(const Phaser::PhaseSet& phasings)
{
    print_phase_sets(std::cout, phasings);
}

} // namespace debug

} // namespace Ocotpus
