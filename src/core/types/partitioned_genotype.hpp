// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef partitioned_genotype_hpp
#define partitioned_genotype_hpp

#include <initializer_list>
#include <array>
#include <utility>
#include <memory>
#include <functional>
#include <iterator>
#include <algorithm>
#include <cassert>
#include <iostream>

#include <boost/functional/hash.hpp>

#include "concepts/equitable.hpp"
#include "utils/append.hpp"
#include "genotype.hpp"

namespace octopus {

template <typename MappableType, std::size_t K>
class PartitionedGenotype
    : public Equitable<PartitionedGenotype<MappableType, K>>, public Mappable<PartitionedGenotype<MappableType, K>>
{
public:
    using MappingDomain = typename Genotype<MappableType>::MappingDomain;
    
    PartitionedGenotype() = default;

    PartitionedGenotype(std::array<Genotype<MappableType>, K> genotype);
    
    PartitionedGenotype(const PartitionedGenotype&)            = default;
    PartitionedGenotype& operator=(const PartitionedGenotype&) = default;
    PartitionedGenotype(PartitionedGenotype&&)                 = default;
    PartitionedGenotype& operator=(PartitionedGenotype&&)      = default;
    
    ~PartitionedGenotype() = default;
    
    const GenomicRegion& mapped_region() const noexcept;

    std::size_t num_partitions() const noexcept;

    const Genotype<MappableType>& partition(std::size_t n) const noexcept;
    Genotype<MappableType>& partition(std::size_t n) noexcept;
    
    const MappableType& operator[](unsigned n) const noexcept;
    MappableType& operator[](unsigned n) noexcept;

    unsigned ploidy() const noexcept;

    friend bool operator==(const PartitionedGenotype<MappableType, K>& lhs, const PartitionedGenotype<MappableType, K>& rhs);
    
private:
    std::array<Genotype<MappableType>, K> partitions_;

    std::pair<std::size_t, unsigned> find_partition(unsigned n) const noexcept;
};

template <typename MappableType, std::size_t K>
PartitionedGenotype<MappableType, K>::PartitionedGenotype(std::array<Genotype<MappableType>, K> genotype)
: partitions_ {std::move(genotype)}
{}

template <typename MappableType, std::size_t K>
const Genotype<MappableType>& PartitionedGenotype<MappableType, K>::partition(std::size_t n) const noexcept
{
    return partitions_[n];
}

template <typename MappableType, std::size_t K>
Genotype<MappableType>& PartitionedGenotype<MappableType, K>::partition(std::size_t n) noexcept
{
    return partitions_[n];
}

template <typename MappableType, std::size_t K>
const MappableType& PartitionedGenotype<MappableType, K>::operator[](unsigned n) const noexcept
{
    const auto p = find_partition(n);
    return partitions_[p.first][p.second];
}

template <typename MappableType, std::size_t K>
MappableType& PartitionedGenotype<MappableType, K>::operator[](unsigned n) noexcept
{
    const auto p = find_partition(n);
    return partitions_[p.first][p.second];
}

template <typename MappableType, std::size_t K>
unsigned PartitionedGenotype<MappableType, K>::ploidy() const noexcept
{
    const static auto add_ploidy = [] (auto total, const auto& genotype) { return total + genotype.ploidy(); };
    return std::accumulate(std::cbegin(partitions_), std::cend(partitions_), 0u, add_ploidy);
}

template <typename MappableType, std::size_t K>
std::pair<std::size_t, unsigned> PartitionedGenotype<MappableType, K>::find_partition(unsigned n) const noexcept
{
    unsigned k {0};
    for (std::size_t i {0}; i < K; ++i) {
        k += partitions_[i].ploidy();
        if (n < k) {
            n = k - n;
            k = i;
            break;
        }
    }
    return {k, n};
}

// free functions

template <typename MappableType, std::size_t K>
bool operator==(const PartitionedGenotype<MappableType, K>& lhs, const PartitionedGenotype<MappableType, K>& rhs)
{
    return lhs.partitions_ == rhs.partitions_;
}

namespace debug {

template <typename S, typename T, std::size_t K>
void print_alleles(S&& stream, const PartitionedGenotype<T, K>& genotype)
{
    for (std::size_t k {0}; k < K - 1; ++k) {
        print_alleles(stream, genotype.partition(k));
        stream << " + ";
    }
    print_alleles(stream, genotype.partition(K - 1));
}

template <typename S, typename T, std::size_t K>
void print_variant_alleles(S&& stream, const PartitionedGenotype<T, K>& genotype)
{
    for (std::size_t k {0}; k < K - 1; ++k) {
        print_variant_alleles(stream, genotype.partition(k));
        stream << " + ";
    }
    print_variant_alleles(stream, genotype.partition(K - 1));
}

template <typename T, std::size_t K>
void print_alleles(const PartitionedGenotype<T, K>& genotype)
{
    print_alleles(std::cout, genotype);
}

template <typename T, std::size_t K>
void print_variant_alleles(const PartitionedGenotype<T, K>& genotype)
{
    print_variant_alleles(std::cout, genotype);
}

} // namespace debug

} // namespace octopus

#endif
