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
#include <iterator>
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
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 0, 4, // N -> A
        4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
    };
    return hash_table[base];
}

template <unsigned char K, typename InputIt>
constexpr auto perfect_kmer_hash(InputIt first)
{
    unsigned k {1};
    
    return std::accumulate(first, std::next(first, K), 0,
                           [&k] (const unsigned curr, const char base) {
                               const auto result = curr + k * perfect_hash(base);
                               k *= 4;
                               return result;
                           });
}

using KmerPerfectHashes = std::vector<std::size_t>;

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
KmerHashTable make_kmer_hash_table(const std::string& sequence)
{
    std::vector<std::vector<std::size_t>> table {num_kmers(K), std::vector<std::size_t> {}};
    
    const auto last_index = sequence.size() - K;
    
    auto it = std::cbegin(sequence);
    
    for (std::size_t index {0}; index <= last_index; ++index, ++it) {
        table[perfect_kmer_hash<K>(it)].push_back(index);
    }
    
    for (auto& bin : table) bin.shrink_to_fit();
    
    return std::make_pair(std::move(table), sequence.size() - K + 1);
}

std::vector<std::size_t>
map_query_to_target(const KmerPerfectHashes& query, const KmerHashTable& target);

template <unsigned char K>
std::vector<std::size_t> map_query_to_target(const std::string& query, const std::string& target)
{
    return map_query_to_target(compute_kmer_hashes<K>(query), make_kmer_hash_table<K>(target));
}

#endif /* kmer_mapping_hpp */
