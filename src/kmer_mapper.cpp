//
//  kmer_mapping.cpp
//  simd_pair_hmm
//
//  Created by Daniel Cooke on 15/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "kmer_mapper.hpp"

void map_query_to_target(const KmerPerfectHashes& query, const KmerHashTable& target,
                         MappedIndexCounts& mapping_counts, std::vector<std::size_t>& result)
{
    map_query_to_target(query, target, mapping_counts, std::back_inserter(result));
}

std::vector<std::size_t>
map_query_to_target(const KmerPerfectHashes& query, const KmerHashTable& target,
                    MappedIndexCounts& mapping_counts)
{
    std::vector<std::size_t> result {};
    result.reserve(1);
    map_query_to_target(query, target, mapping_counts, result);
    return result;
}

std::vector<std::size_t>
map_query_to_target(const KmerPerfectHashes& query, const KmerHashTable& target)
{
    MappedIndexCounts mapping_counts(target.second, 0);
    return  map_query_to_target(query, target, mapping_counts);
}
