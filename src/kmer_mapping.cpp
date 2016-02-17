//
//  kmer_mapping.cpp
//  simd_pair_hmm
//
//  Created by Daniel Cooke on 15/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "kmer_mapping.hpp"

std::vector<std::size_t>
extract_maximum_hash_hit_indicies(const KmerPerfectHashes& query, const KmerHashTable& target,
                                  const std::size_t target_size)
{
    std::vector<unsigned> hit_counts(target_size, 0);
    
    unsigned max_hit_count {0};
    
    for (std::size_t query_index {0}; query_index < query.size(); ++query_index) {
        for (const auto target_index : target[query[query_index]]) {
            if (target_index >= query_index) {
                const auto mapping_begin = target_index - query_index;
                
                if (++hit_counts[mapping_begin] > max_hit_count) {
                    max_hit_count = hit_counts[mapping_begin];
                }
            }
        }
    }
    
    std::vector<std::size_t> result {};
    
    if (max_hit_count > 0) {
        result.reserve(1);
        
        for (std::size_t i {0}; i < target_size; ++i) {
            if (hit_counts[i] == max_hit_count) {
                result.push_back(i);
            }
        }
        
        result.shrink_to_fit();
    }
    
    return result;
}
