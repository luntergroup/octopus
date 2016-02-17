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
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
    };
    return hash_table[base];
}

template <unsigned char K, typename InputIt>
constexpr auto perfect_kmer_hash(InputIt first)
{
    unsigned k {1}, result {0};
    
    const auto last = std::next(first, K);
    
    for (; first != last; ++first) {
        result += k * perfect_hash(*first);
        k *= 4;
    }
    
    return result;
}

using KmerPerfectHashes = std::vector<std::size_t>;

template <unsigned char K>
auto compute_kmer_hashes(const std::string& sequence)
{
    KmerPerfectHashes result(sequence.size() - K);
    
    auto result_it = std::begin(result);
    
    auto it = std::cbegin(sequence);
    
    const auto last = std::prev(std::cend(sequence), K);
    
    for (; it != last; ++it, ++result_it) {
        *result_it = perfect_kmer_hash<K>(it);
    }
    
    return result;
}

using KmerHashTable = std::vector<std::vector<std::size_t>>;

template <unsigned char K>
auto make_kmer_hash_table(const std::string& sequence)
{
    KmerHashTable result {num_kmers(K), KmerHashTable::value_type {}};
    
    const auto last_index = sequence.size() - K;
    
    auto it = std::cbegin(sequence);
    
    for (std::size_t index {0}; index < last_index; ++index, ++it) {
        result[perfect_kmer_hash<K>(it)].push_back(index);
    }
    
    for (auto& bin : result) bin.shrink_to_fit();
    
    return result;
}

std::vector<std::size_t>
extract_maximum_hash_hit_indicies(const KmerPerfectHashes& query, const KmerHashTable& target,
                                  std::size_t target_size);

template <unsigned char K>
std::vector<std::size_t> map_query_to_target(const std::string& query, const std::string& target)
{
    const auto query_hashes      = compute_kmer_hashes<K>(query);
    const auto target_hash_table = make_kmer_hash_table<K>(target);    
    return extract_maximum_hash_hit_indicies(query_hashes, target_hash_table, target.size() - K);
}

#endif /* kmer_mapping_hpp */
