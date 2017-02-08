// Copyright (c) 2016 Daniel Cooke
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

Phaser::Phaser(Phred<double> min_phase_score) : min_phase_score_ {min_phase_score} {}

namespace {

using GenotypeReference  = std::reference_wrapper<const Genotype<Haplotype>>;
using PhaseComplementSet = std::deque<GenotypeReference>;
using PartitionIterator  = std::vector<GenomicRegion>::const_iterator;

struct PhaseComplementHash
{
    PhaseComplementHash(PartitionIterator first, PartitionIterator last) : first_ {first}, last_ {last} {}
    
    auto operator()(const GenotypeReference& genotype) const
    {
        using boost::hash_combine;
        std::size_t result {0};
        std::for_each(first_, last_, [&] (const auto& region) {
            hash_combine(result, std::hash<Genotype<Haplotype>>()(splice<Haplotype>(genotype.get(), region)));
        });
        return result;
    }
private:
    PartitionIterator first_, last_;
};

struct PhaseComplementEqual
{
    PhaseComplementEqual(PartitionIterator first, PartitionIterator last)
    : first_ {first}
    , last_ {last}
    {}
    
    bool operator()(const GenotypeReference& lhs, const GenotypeReference& rhs) const
    {
        return std::all_of(first_, last_,
                           [&] (const auto& region) {
                               return are_equal_in_region<Haplotype>(lhs.get(), rhs.get(), region);
                           });
    }
private:
    PartitionIterator first_, last_;
};

using PhaseComplementSetMap = std::unordered_map<GenotypeReference, PhaseComplementSet, PhaseComplementHash, PhaseComplementEqual>;
using PhaseComplementSets = std::vector<PhaseComplementSet>;

template <typename Container>
auto make_phase_completement_set_map(const Container& genotypes,
                                     PartitionIterator first, PartitionIterator last)
{
    PhaseComplementSetMap result {genotypes.size(), PhaseComplementHash {first, last},
                                    PhaseComplementEqual {first, last}};
    result.reserve(genotypes.size());
    return result;
}

void insert(const Genotype<Haplotype>& genotype, PhaseComplementSetMap& curr_result)
{
    const auto itr = curr_result.find(genotype);
    if (itr == std::end(curr_result)) {
        curr_result.emplace(std::piecewise_construct,
                            std::forward_as_tuple(genotype),
                            std::forward_as_tuple(1, genotype));
    } else {
        itr->second.emplace_back(genotype);
    }
}

template <typename Container>
PhaseComplementSets
generate_phase_complement_sets(const Container& genotypes,
                               PartitionIterator first, PartitionIterator last)
{
    auto complement_sets = make_phase_completement_set_map(genotypes, first, last);
    complement_sets.emplace(std::piecewise_construct,
                            std::forward_as_tuple(genotypes.front()),
                            std::forward_as_tuple(1, genotypes.front()));
    std::for_each(std::next(std::cbegin(genotypes)), std::cend(genotypes),
                  [&] (const auto& genotype) { insert(genotype, complement_sets); });
    PhaseComplementSets result {};
    result.reserve(complement_sets.size());
    for (auto&& p : complement_sets) {
        result.emplace_back(std::move(p.second));
    }
    return result;
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
    return std::max(0.0, -std::accumulate(std::cbegin(phase_set), std::cend(phase_set), 0.0,
                                          [&genotype_posteriors, norm] (const auto curr, const auto& genotype) {
                                              const auto p = genotype_posteriors.at(genotype) / norm;
                                              return curr + p * std::log2(p);
                                          }));
}

double maximum_entropy(const size_t num_elements)
{
    return std::log2(num_elements);
}

template <typename Map>
double calculate_relative_entropy(const PhaseComplementSet& phase_set, const Map& genotype_posteriors)
{
    if (phase_set.size() < 2) return 1.0;
    return 1.0 - calculate_entropy(phase_set, genotype_posteriors) / maximum_entropy(phase_set.size());
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

auto min_phase_score(const double p, const std::size_t num_genotypes)
{
    if (maths::almost_one(p) || num_genotypes == 1) {
        static const Phred<double> max_possible_score {Phred<double>::Probability {0.0}};
        return max_possible_score;
    } else {
        const auto s = p * std::log2(p) + (1.0 - p) * std::log2((1.0 - p) / (num_genotypes - 1));
        const auto y = 1.0 + s / maximum_entropy(2);
        return Phred<double> {Phred<double>::Probability {y}};
    }
}

auto min_phase_score(const Genotype<Haplotype>& called_genotype, const Phaser::SampleGenotypePosteriorMap& genotype_posteriors)
{
    return min_phase_score(genotype_posteriors[called_genotype], genotype_posteriors.size());
}

std::vector<GenotypeReference> extract_genotypes(const Phaser::GenotypePosteriorMap& genotype_posteriors)
{
    return extract_key_refs(genotype_posteriors);
}

using GenotypeSplicePosteriorMap = std::unordered_map<Genotype<Haplotype>, double>;

template <typename Container>
auto splice_and_marginalise(const Container& genotypes,
                            const Phaser::SampleGenotypePosteriorMap& genotype_posteriors,
                            const GenomicRegion& region)
{
    GenotypeSplicePosteriorMap splice_posteriors {genotypes.size()};
    std::vector<Genotype<Haplotype>> splices {};
    splices.reserve(genotypes.size());
    for (const auto& p : genotype_posteriors) {
        auto splice_genotype = splice<Haplotype>(p.first, region);
        splice_posteriors[splice_genotype] += p.second;
        splices.push_back(std::move(splice_genotype));
    }
    std::sort(std::begin(splices), std::end(splices), GenotypeLess {});
    splices.erase(std::unique(std::begin(splices), std::end(splices)), std::end(splices));
    return std::make_pair(std::move(splices), std::move(splice_posteriors));
}

} // namespace

boost::optional<Phaser::PhaseSet>
Phaser::try_phase(const std::vector<Haplotype>& haplotypes,
                  const GenotypePosteriorMap& genotype_posteriors,
                  const std::vector<GenomicRegion>& regions) const
{
    // TODO
    return boost::none;
}

Phaser::PhaseSet::SamplePhaseRegions
force_phase_sample(const GenomicRegion& region,
                   const std::vector<GenomicRegion>& partitions,
                   const std::vector<GenotypeReference>& genotypes,
                   const Phaser::SampleGenotypePosteriorMap& genotype_posteriors,
                   const Phred<double> min_phase_score)
{
    auto first_partition = std::cbegin(partitions);
    auto last_partition  = std::cend(partitions);
    auto phase_set = generate_phase_complement_sets(genotypes, first_partition, last_partition);
    auto phase_score = calculate_phase_score(phase_set, genotype_posteriors);
    if (phase_score >= min_phase_score) {
        return {Phaser::PhaseSet::PhaseRegion {region, phase_score}};
    }
    Phaser::PhaseSet::SamplePhaseRegions result {};
    --last_partition;
    std::vector<Genotype<Haplotype>> splices;
    GenotypeSplicePosteriorMap splice_posteriors;
    while (first_partition != std::cend(partitions)) {
        auto curr_region = encompassing_region(first_partition, last_partition);
        std::tie(splices, splice_posteriors) = splice_and_marginalise(genotypes, genotype_posteriors, curr_region);
        phase_set = generate_phase_complement_sets(splices, first_partition, last_partition);
        phase_score = calculate_phase_score(phase_set, splice_posteriors);
        if (phase_score >= min_phase_score || std::distance(first_partition, last_partition) == 1) {
            result.emplace_back(encompassing_region(first_partition, last_partition), phase_score);
            first_partition = last_partition;
            last_partition  = std::cend(partitions);
        } else {
            --last_partition;
        }
    }
    return result;
}

Phaser::PhaseSet
Phaser::force_phase(const std::vector<Haplotype>& haplotypes,
                    const GenotypePosteriorMap& genotype_posteriors,
                    const std::vector<GenomicRegion>& regions,
                    boost::optional<GenotypeCallMap> genotype_calls) const
{
    assert(!haplotypes.empty());
    assert(!genotype_posteriors.empty1());
    assert(!genotype_posteriors.empty2());
    assert(std::is_sorted(std::cbegin(regions), std::cend(regions)));
    const auto& haplotype_region = haplotypes.front().mapped_region();
    const auto genotypes = extract_genotypes(genotype_posteriors);
    const auto partitions = extract_covered_regions(regions);
    PhaseSet result {haplotype_region};
    result.phase_regions.reserve(genotype_posteriors.size1());
    if (genotypes.front().get().ploidy() == 1 || partitions.size() == 1) {
        for (const auto& p : genotype_posteriors) {
            if (max_phase_score_) {
                result.phase_regions[p.first].emplace_back(haplotype_region, *max_phase_score_);
            } else {
                static const Phred<double> max_possible_score {Phred<double>::Probability {0.0}};
                result.phase_regions[p.first].emplace_back(haplotype_region, max_possible_score);
            }
        }
        return result;
    }
    for (const auto& p : genotype_posteriors) {
        if (genotype_calls && max_phase_score_ && min_phase_score(genotype_calls->at(p.first), p.second) >= *max_phase_score_) {
            result.phase_regions[p.first].emplace_back(haplotype_region, *max_phase_score_);
        } else {
            auto phases = force_phase_sample(haplotype_region, partitions, genotypes, p.second, min_phase_score_);
            if (max_phase_score_) {
                for (auto& phase : phases) phase.score = std::min(phase.score, *max_phase_score_);
            }
            result.phase_regions.emplace(p.first, std::move(phases));
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
