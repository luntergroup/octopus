//
//  kmer_mapper.cpp
//  Octopus
//
//  Created by Daniel Cooke on 23/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "kmer_mapper.hpp"

#include "maths.hpp"

namespace Octopus
{
    KmerMapper::KmerMapper(const ReadMap& reads, const std::vector<Haplotype>& haplotypes)
    :
    read_cache_ {Maths::sum_sizes(reads)},
    haplotype_cache_ {haplotypes.size()}
    {
        for (const auto& s : reads) {
            for (const auto& read : s.second) {
                if (read_cache_.count(read.get_sequence()) == 0) {
                    read_cache_.emplace(read.get_sequence(),
                                        compute_kmer_hashes<KMER_SIZE>(read.get_sequence()));
                }
            }
        }
        
        for (const auto& haplotype : haplotypes) {
            if (haplotype_cache_.count(haplotype.get_sequence()) == 0) {
                haplotype_cache_.emplace(haplotype.get_sequence(),
                                         make_kmer_hash_table<KMER_SIZE>(haplotype.get_sequence()));
            }
        }
    }
    
    std::vector<std::size_t> KmerMapper::map(const AlignedRead& read, const Haplotype& haplotype) const
    {
        return map_query_to_target(read_cache_.at(read.get_sequence()),
                                   haplotype_cache_.at(haplotype.get_sequence()));
    }
    
    void KmerMapper::clear() noexcept
    {
        read_cache_.clear();
        haplotype_cache_.clear();
    }
} // namespace Octopus
