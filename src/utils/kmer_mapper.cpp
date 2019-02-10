// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "kmer_mapper.hpp"

namespace octopus {

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

} // namespace octopus
