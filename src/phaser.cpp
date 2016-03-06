//
//  phaser.cpp
//  Octopus
//
//  Created by Daniel Cooke on 20/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "phaser.hpp"

#include <deque>
#include <algorithm>
#include <numeric>
#include <cstddef>
#include <cmath>
#include <functional>
#include <utility>
#include <cassert>

#include <boost/functional/hash.hpp>

#include "mappable.hpp"
#include "mappable_algorithms.hpp"

#include <iostream> // DEBUG

namespace Octopus
{
    Phaser::Phaser(const double min_phase_score)
    :
    min_phase_score_ {min_phase_score}
    {}
    
    namespace
    {
        using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;
        
        using PhaseComplementSet = std::deque<GenotypeReference>;
        
        using PartitionIterator = std::vector<GenomicRegion>::const_iterator;
        
        struct PhaseComplementHash
        {
            explicit PhaseComplementHash(PartitionIterator first, PartitionIterator last)
            : first_ {first}, last_ {last} {}
            
            std::size_t operator()(const GenotypeReference& genotype) const
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
            explicit PhaseComplementEqual(PartitionIterator first, PartitionIterator last)
            : first_ {first}, last_ {last} {}
            
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
        
        using PhaseComplementSetMap = std::unordered_map<GenotypeReference, PhaseComplementSet,
                                                            PhaseComplementHash,
                                                            PhaseComplementEqual>;
        
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
        
        void add_to_phase_complement_set(const Genotype<Haplotype>& genotype,
                                         PhaseComplementSetMap& curr_result)
        {
            const auto it = curr_result.find(genotype);
            if (it == std::end(curr_result)) {
                curr_result.emplace(std::piecewise_construct,
                                    std::forward_as_tuple(genotype),
                                    std::forward_as_tuple(1, genotype));
            } else {
                it->second.emplace_back(genotype);
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
                          [&] (const auto& genotype) {
                              add_to_phase_complement_set(genotype, complement_sets);
                          });
            
            PhaseComplementSets result {};
            result.reserve(complement_sets.size());
            
            for (auto&& p : complement_sets) {
                result.emplace_back(std::move(p.second));
            }
            
            return result;
        }
        
        template <typename Map>
        double marginalise(const PhaseComplementSet& phase_set,
                           const Map& genotype_posteriors)
        {
            return std::accumulate(std::cbegin(phase_set), std::cend(phase_set), 0.0,
                                   [&] (const auto curr, const auto& genotype) {
                                       return curr + genotype_posteriors.at(genotype);
                                   });
        }
        
        template <typename Map>
        double calculate_entropy(const PhaseComplementSet& phase_set,
                                 const Map& genotype_posteriors)
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
        double calculate_relative_entropy(const PhaseComplementSet& phase_set,
                                          const Map& genotype_posteriors)
        {
            if (phase_set.size() < 2) return 1.0;
            return 1.0 - calculate_entropy(phase_set, genotype_posteriors) / maximum_entropy(phase_set.size());
        }
        
        template <typename Map>
        double calculate_phase_score(const PhaseComplementSet& phase_set,
                                     const Map& genotype_posteriors)
        {
            return marginalise(phase_set, genotype_posteriors) * calculate_relative_entropy(phase_set, genotype_posteriors);
        }
        
        template <typename Map>
        double calculate_phase_score(const PhaseComplementSets& phase_sets,
                                     const Map& genotype_posteriors)
        {
            return std::accumulate(std::cbegin(phase_sets), std::cend(phase_sets), 0.0,
                                   [&] (const auto curr, const auto& phase_set) {
                                       return curr + calculate_phase_score(phase_set, genotype_posteriors);
                                   });
        }
        
        std::vector<GenotypeReference>
        extract_genotypes(const Phaser::GenotypePosteriorMap& genotype_posteriors)
        {
            return extract_key_refs(genotype_posteriors);
        }
        
        using GenotypeSplicePosteriorMap = std::unordered_map<Genotype<Haplotype>, double>;
        
        
        
        template <typename Container>
        auto splice_and_marginalise(const Container& genotypes,
                                    const Phaser::SampleGenotypePosteriorMap& genotype_posteriors,
                                    const GenomicRegion& region)
        {
            auto splices = splice_all<Haplotype>(genotypes, region);
            
            GenotypeSplicePosteriorMap splice_posteriors {splices.size()};
            
            for (const auto& splice : splices) {
                splice_posteriors.emplace(splice, 0.0);
            }
            
            for (const auto& p : genotype_posteriors) {
                splice_posteriors.at(splice<Haplotype>(p.first, region)) += p.second;
            }
            
            return std::make_pair(std::move(splices), std::move(splice_posteriors));
        }
    } // namespace
    
    boost::optional<Phaser::PhaseSet>
    Phaser::try_phase(const std::vector<Haplotype>& haplotypes,
                      const GenotypePosteriorMap& genotype_posteriors,
                      const std::vector<Variant>& candidates)
    {
        return boost::none;
    }
    
    Phaser::PhaseSet::SamplePhaseRegions
    force_phase_sample(const GenomicRegion& region, const std::vector<GenomicRegion>& partitions,
                       const std::vector<GenotypeReference>& genotypes,
                       const Phaser::SampleGenotypePosteriorMap& genotype_posteriors,
                       const double min_phase_score)
    {
        auto first_partition = std::cbegin(partitions);
        auto last_partition  = std::cend(partitions);
        
        auto phase_set = generate_phase_complement_sets(genotypes, first_partition, last_partition);
        
        auto phase_score = calculate_phase_score(phase_set, genotype_posteriors);
        
        if (phase_score >= min_phase_score) {
            return {Phaser::PhaseSet::PhaseRegion {region, phase_score}};
        }
        
        Phaser::PhaseSet::SamplePhaseRegions result {};
        
        last_partition = std::next(first_partition);
        
        double previous_phase_score {0};
        
        while (first_partition != std::cend(partitions)) {
            auto curr_region = encompassing_region(first_partition, last_partition);
            
            auto splice_posteriors = splice_and_marginalise(genotypes, genotype_posteriors, curr_region);
            
            phase_set = generate_phase_complement_sets(splice_posteriors.first,
                                                       first_partition, last_partition);
            
            phase_score = calculate_phase_score(phase_set, splice_posteriors.second);
            
            if (phase_score + 0.1 < previous_phase_score) {
                auto phase_region = encompassing_region(first_partition, std::prev(last_partition));
                result.emplace_back(std::move(phase_region), previous_phase_score);
                first_partition = std::prev(last_partition);
                previous_phase_score = 0.0;
            } else {
                if (last_partition == std::cend(partitions)) {
                    result.emplace_back(curr_region, phase_score);
                    first_partition = last_partition;
                } else {
                    previous_phase_score = phase_score;
                    ++last_partition;
                }
            }
        }
        
        return result;
    }
    
    Phaser::PhaseSet
    Phaser::force_phase(const std::vector<Haplotype>& haplotypes,
                        const GenotypePosteriorMap& genotype_posteriors,
                        const std::vector<Variant>& candidates)
    {
        assert(!haplotypes.empty());
        assert(!candidates.empty());
        assert(!genotype_posteriors.empty1());
        assert(!genotype_posteriors.empty2());
        
        const auto& haplotype_region = haplotypes.front().get_region();
        
        PhaseSet result {haplotype_region};
        
        result.phase_regions.reserve(genotype_posteriors.size1());
        
        const auto genotypes = extract_genotypes(genotype_posteriors);
        
        for (const auto& p : genotype_posteriors) {
            result.phase_regions[p.first].emplace_back(haplotype_region, 1);
        }
        return result;
        
//        if (genotypes.front().get().ploidy() == 1 || candidates.size() == 1) {
//            for (const auto& p : genotype_posteriors) {
//                result.phase_regions[p.first].emplace_back(haplotype_region, 1);
//            }
//            return result;
//        }
//        
//        const auto partitions = extract_covered_regions(candidates);
//        
//        for (const auto& p : genotype_posteriors) {
//            result.phase_regions.emplace(p.first, force_phase_sample(haplotype_region, partitions,
//                                                                     genotypes, p.second,
//                                                                     min_phase_score_));
//        }
//        
//        return result;
    }
    
} // namespace Ocotpus
