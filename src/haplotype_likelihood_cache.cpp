//
//  haplotype_likelihood_cache.cpp
//  Octopus
//
//  Created by Daniel Cooke on 25/11/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "haplotype_likelihood_cache.hpp"

#include <vector>
#include <algorithm>
#include <iostream>
#include <utility>

#include "kmer_mapper.hpp"

#include "timers.hpp"

namespace Octopus
{
// public methods

HaplotypeLikelihoodCache::HaplotypeLikelihoodCache(const ReadSet& reads,
                                                   const std::vector<Haplotype>& haplotypes,
                                                   HaplotypeLikelihoodModel::InactiveRegionState flank_state)
:
HaplotypeLikelihoodCache {HaplotypeLikelihoodModel {KmerMapper(reads, haplotypes)}, reads, haplotypes, flank_state}
{}

HaplotypeLikelihoodCache::HaplotypeLikelihoodCache(HaplotypeLikelihoodModel error_model,
                                                   const ReadSet& reads,
                                                   const std::vector<Haplotype>& haplotypes,
                                                   HaplotypeLikelihoodModel::InactiveRegionState flank_state)
:
error_model_ {std::move(error_model)}
{    
    cache_.assign_keys(std::cbegin(reads), std::cend(reads));
    
    cache_.reserve1(haplotypes.size());
    
    std::vector<double> probabilities(reads.size());
    
    for (const auto& haplotype : haplotypes) {
        std::transform(std::cbegin(reads), std::cend(reads), std::begin(probabilities),
                       [this, &haplotype, flank_state] (const auto& read) {
                           return error_model_.log_probability(read, haplotype, flank_state);
                       });
        
        cache_.insert_at(haplotype, std::cbegin(probabilities), std::cend(probabilities));
    }
}

double HaplotypeLikelihoodCache::log_probability(const AlignedRead& read, const Haplotype& haplotype) const
{
    resume_timer(likelihood_cache_timer);
    const auto result = cache_(haplotype, read);
    pause_timer(likelihood_cache_timer);
    return result;
}

void HaplotypeLikelihoodCache::clear()
{
    cache_.clear();
}

namespace debug
{
    void print_read_haplotype_liklihoods(const std::vector<Haplotype>& haplotypes,
                                         const ReadMap& reads,
                                         HaplotypeLikelihoodCache& haplotype_likelihoods,
                                         size_t n)
    {
        auto m = std::min(n, haplotypes.size());
        
        std::cout << "debug: top " << m << " haplotype likelihoods for each read in each sample" << std::endl;
        
        using HaplotypeReference = std::reference_wrapper<const Haplotype>;
        
        for (const auto& sample_reads : reads) {
            std::cout << "Sample: " << sample_reads.first << ":" << std::endl;
            
            for (const auto& read : sample_reads.second) {
                std::cout << "\tRead: " << read.get_region() << " " << read.get_cigar_string() << ":" << std::endl;
                
                std::vector<std::pair<HaplotypeReference, double>> top {};
                top.reserve(haplotypes.size());
                
                for (const auto& haplotype : haplotypes) {
                    top.emplace_back(haplotype, haplotype_likelihoods.log_probability(read, haplotype));
                }
                
                std::sort(std::begin(top), std::end(top),
                          [] (const auto& lhs, const auto& rhs) {
                              return lhs.second > rhs.second;
                          });
                
                for (unsigned i {}; i < m; ++i) {
                    std::cout << "\t\t ";
                    print_variant_alleles(top[i].first);
                    std::cout << " " << std::setprecision(10) << top[i].second << std::endl;
                }
            }
        }
    }
} // namespace debug
} // namespace Octopus
