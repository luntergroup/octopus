//
//  kmer_mapping.hpp
//  simd_pair_hmm
//
//  Created by Daniel Cooke on 15/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef kmer_mapping_hpp
#define kmer_mapping_hpp

#include <string>
#include <vector>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <algorithm>
#include <numeric>

constexpr auto num_kmers(const unsigned char k) noexcept
{
    return 2 << (2 * static_cast<std::size_t>(k) - 1);
}

template <typename T = short>
constexpr auto perfect_hash(const char base) noexcept
{
    constexpr std::array<T, 128> hash_table
    {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    return hash_table[base];
}

using KmerHashType = std::uint_fast32_t;

template <unsigned char K, typename InputIt>
constexpr auto perfect_kmer_hash(InputIt first)
{
    unsigned k {1};
    
    return std::accumulate(first, std::next(first, K), KmerHashType {0},
                           [&k] (const unsigned curr, const char base) {
                               const auto result = curr + k * perfect_hash(base);
                               k *= 4;
                               return result;
                           });
}

using KmerPerfectHashes = std::vector<KmerHashType>;

template <unsigned char K>
auto compute_kmer_hashes(const std::string& sequence)
{
    KmerPerfectHashes result(sequence.size() - K + 1);
    
    auto result_it = std::begin(result);
    
    for (auto it = std::cbegin(sequence); it != std::prev(std::cend(sequence), K - 1); ++it, ++result_it) {
        *result_it = perfect_kmer_hash<K>(it);
    }
    
    return result;
}

using KmerHashTable = std::pair<std::vector<std::vector<std::size_t>>, std::size_t>;

template <unsigned char K>
KmerHashTable init_kmer_hash_table()
{
    return std::make_pair(std::vector<std::vector<std::size_t>> {num_kmers(K), std::vector<std::size_t> {}}, 0);
}

inline void clear_kmer_hash_table(KmerHashTable& table)
{
    for (auto& bin : table.first) bin.clear();
    table.second = 0;
}

template <unsigned char K>
void populate_kmer_hash_table(const std::string& sequence, KmerHashTable& result)
{
    const auto last_index = sequence.size() - K;
    
    auto it = std::cbegin(sequence);
    
    for (std::size_t index {0}; index <= last_index; ++index, ++it) {
        result.first[perfect_kmer_hash<K>(it)].push_back(index);
    }
    
    for (auto& bin : result.first) bin.shrink_to_fit();
    
    result.second = sequence.size() - K + 1;
}

template <unsigned char K>
KmerHashTable make_kmer_hash_table(const std::string& sequence)
{
    auto result = init_kmer_hash_table<K>();
    populate_kmer_hash_table<K>(sequence, result);
    return result;
}

using MappedIndexCounts = std::vector<unsigned>;

inline MappedIndexCounts init_mapping_counts(const KmerHashTable& target)
{
    return MappedIndexCounts(target.second, 0);
}

inline void reset_mapping_counts(MappedIndexCounts& mapping_counts)
{
    std::fill(std::begin(mapping_counts), std::end(mapping_counts), 0);
}

template <typename OutputIt>
OutputIt map_query_to_target(const KmerPerfectHashes& query, const KmerHashTable& target,
                             MappedIndexCounts& mapping_counts, OutputIt result)
{
    unsigned max_hit_count {0};
    std::size_t first_max_hit_index {};
    unsigned num_max_hits {0};
    
    for (std::size_t query_index {0}; query_index < query.size(); ++query_index) {
        for (const auto target_index : target.first[query[query_index]]) {
            if (target_index >= query_index) {
                const auto mapping_begin = target_index - query_index;
                
                if (++mapping_counts[mapping_begin] > max_hit_count) {
                    max_hit_count = mapping_counts[mapping_begin];
                    first_max_hit_index = mapping_begin;
                    num_max_hits = 1;
                } else if (mapping_counts[mapping_begin] == max_hit_count) {
                    ++num_max_hits;
                    
                    if (mapping_begin < first_max_hit_index) {
                        first_max_hit_index = mapping_begin;
                    }
                }
            }
        }
    }
    
    if (max_hit_count > 0) {
        *result++ = first_max_hit_index++;
        
        --num_max_hits;
        
        while (num_max_hits > 0) {
            if (mapping_counts[first_max_hit_index] == max_hit_count) {
                *result++ = first_max_hit_index;
                --num_max_hits;
            }
            ++first_max_hit_index;
        }
    }
    
    return result;
}

void map_query_to_target(const KmerPerfectHashes& query, const KmerHashTable& target,
                         MappedIndexCounts& mapping_counts, std::vector<std::size_t>& result);

std::vector<std::size_t>
map_query_to_target(const KmerPerfectHashes& query, const KmerHashTable& target,
                    MappedIndexCounts& mapping_counts);

std::vector<std::size_t>
map_query_to_target(const KmerPerfectHashes& query, const KmerHashTable& target);

template <unsigned char K>
std::vector<std::size_t> map_query_to_target(const std::string& query, const std::string& target)
{
    return map_query_to_target(compute_kmer_hashes<K>(query), make_kmer_hash_table<K>(target));
}

#endif /* kmer_mapping_hpp */
