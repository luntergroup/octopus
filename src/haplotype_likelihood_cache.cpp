//
//  haplotype_likelihood_cache.cpp
//  Octopus
//
//  Created by Daniel Cooke on 25/11/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "haplotype_likelihood_cache.hpp"

#include <algorithm>
#include <iostream>
#include <vector>
#include <utility>

#include "maths.hpp"

namespace Octopus
{
// public methods

HaplotypeLikelihoodCache::HaplotypeLikelihoodCache(const ReadMap& reads,
                                                   const std::vector<Haplotype>& haplotypes)
:
error_model_ {},
cache_ {Maths::sum_sizes(reads)},
max_num_haplotypes_ {haplotypes.size()}
{
    for (const auto& sample_reads : reads) {
        for (const auto& read : sample_reads.second) {
            for (const auto& haplotype : haplotypes) {
                cache(read, haplotype, error_model_.log_probability(read, haplotype));
            }
        }
    }
}

HaplotypeLikelihoodCache::HaplotypeLikelihoodCache(HaplotypeLikelihoodModel error_model,
                                                   const ReadMap& reads,
                                                   const std::vector<Haplotype>& haplotypes)
:
error_model_ {std::move(error_model)},
cache_ {Maths::sum_sizes(reads)},
max_num_haplotypes_ {haplotypes.size()}
{
    for (const auto& sample_reads : reads) {
        for (const auto& read : sample_reads.second) {
            for (const auto& haplotype : haplotypes) {
                cache(read, haplotype, error_model_.log_probability(read, haplotype));
            }
        }
    }
}

double HaplotypeLikelihoodCache::log_probability(const AlignedRead& read, const Haplotype& haplotype) const
{
    return get_cached(read, haplotype);
}

void HaplotypeLikelihoodCache::clear()
{
    cache_.clear();
}

// private methods

bool HaplotypeLikelihoodCache::is_cached(const AlignedRead& read, const Haplotype& haplotype) const noexcept
{
    return cache_.count(read) == 1 && cache_.at(read).count(haplotype) == 1;
}

void HaplotypeLikelihoodCache::cache(const AlignedRead& read, const Haplotype& haplotype, double value) const
{
    if (cache_.count(read) == 0) {
        cache_.emplace(std::piecewise_construct, std::forward_as_tuple(read),
                       std::forward_as_tuple(max_num_haplotypes_));
    }
    cache_[read].emplace(haplotype, value);
}

double HaplotypeLikelihoodCache::get_cached(const AlignedRead& read, const Haplotype& haplotype) const
{
    return cache_.at(read).at(haplotype);
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
                    std::cout << "\t\t* ";
                    print_variant_alleles(top[i].first);
                    std::cout << " " << std::setprecision(10) << top[i].second << std::endl;
                }
            }
        }
    }
} // namespace debug
} // namespace Octopus
