//
//  haplotype_likelihood_cache.cpp
//  Octopus
//
//  Created by Daniel Cooke on 25/11/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "haplotype_likelihood_cache.hpp"

#include "maths.hpp"

namespace Octopus
{
    // public methods
    
    HaplotypeLikelihoodCache::HaplotypeLikelihoodCache(const ReadMap& reads,
                                                       const std::vector<Haplotype>& haplotypes)
    :
    read_model_ {},
    cache_ {},
    max_num_reads_ {Maths::sum_sizes(reads)},
    max_num_haplotypes_ {haplotypes.size()}
    {
        cache_.reserve(max_num_reads_);
        for (const auto& sample_reads : reads) {
            for (const auto& read : sample_reads.second) {
                cache_[read].reserve(max_num_haplotypes_);
                for (const auto& haplotype : haplotypes) {
                    cache(read, haplotype, read_model_.log_probability(read, haplotype));
                }
            }
        }
    }
    
    HaplotypeLikelihoodCache::HaplotypeLikelihoodCache(SingleReadModel read_model,
                                                       const ReadMap& reads,
                                                       const std::vector<Haplotype>& haplotypes)
    :
    read_model_ {std::move(read_model)},
    cache_ {},
    max_num_reads_ {Maths::sum_sizes(reads)},
    max_num_haplotypes_ {haplotypes.size()}
    {
        cache_.reserve(max_num_reads_);
        for (const auto& sample_reads : reads) {
            for (const auto& read : sample_reads.second) {
                cache_[read].reserve(max_num_haplotypes_);
                for (const auto& haplotype : haplotypes) {
                    cache(read, haplotype, read_model_.log_probability(read, haplotype));
                }
            }
        }
    }
    
    double HaplotypeLikelihoodCache::log_probability(const AlignedRead& read, const Haplotype& haplotype)
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
    
    void HaplotypeLikelihoodCache::cache(const AlignedRead& read, const Haplotype& haplotype, double value)
    {
        if (cache_[read].empty()) {
            cache_[read].reserve(max_num_haplotypes_);
        }
        cache_[read].emplace(haplotype, value);
    }
    
    double HaplotypeLikelihoodCache::get_cached(const AlignedRead& read, const Haplotype& haplotype) const
    {
        return cache_.at(read).at(haplotype);
    }
} // namespace Octopus
