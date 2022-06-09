// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef partitioned_genotype_hpp
#define partitioned_genotype_hpp

#include <initializer_list>
#include <vector>
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

template <typename MappableType>
class PartitionedGenotype
    : public Equitable<PartitionedGenotype<MappableType>>, public Mappable<PartitionedGenotype<MappableType>>
{
public:
    using MappingDomain = typename Genotype<MappableType>::MappingDomain;
    
    PartitionedGenotype() = default;

    PartitionedGenotype(std::size_t k);
    PartitionedGenotype(std::vector<Genotype<MappableType>> genotype);
    
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

    friend bool operator==(const PartitionedGenotype& lhs, const PartitionedGenotype& rhs)
    {
        return lhs.partitions_ == rhs.partitions_;
    }
    
private:
    std::vector<Genotype<MappableType>> partitions_;

    std::pair<std::size_t, unsigned> find_partition(unsigned n) const noexcept;
};

template <typename MappableType>
PartitionedGenotype<MappableType>::PartitionedGenotype(std::size_t n)
: partitions_(n)
{}

template <typename MappableType>
PartitionedGenotype<MappableType>::PartitionedGenotype(std::vector<Genotype<MappableType>> genotype)
: partitions_ {std::move(genotype)}
{}

template <typename MappableType>
std::size_t PartitionedGenotype<MappableType>::num_partitions() const noexcept
{
    return partitions_.size();
}

template <typename MappableType>
const Genotype<MappableType>& PartitionedGenotype<MappableType>::partition(std::size_t n) const noexcept
{
    return partitions_[n];
}

template <typename MappableType>
Genotype<MappableType>& PartitionedGenotype<MappableType>::partition(std::size_t n) noexcept
{
    return partitions_[n];
}

template <typename MappableType>
const MappableType& PartitionedGenotype<MappableType>::operator[](unsigned n) const noexcept
{
    const auto p = find_partition(n);
    return partitions_[p.first][p.second];
}

template <typename MappableType>
MappableType& PartitionedGenotype<MappableType>::operator[](unsigned n) noexcept
{
    const auto p = find_partition(n);
    return partitions_[p.first][p.second];
}

template <typename MappableType>
unsigned PartitionedGenotype<MappableType>::ploidy() const noexcept
{
    const static auto add_ploidy = [] (auto total, const auto& genotype) { return total + genotype.ploidy(); };
    return std::accumulate(std::cbegin(partitions_), std::cend(partitions_), 0u, add_ploidy);
}

template <typename MappableType>
std::pair<std::size_t, unsigned> PartitionedGenotype<MappableType>::find_partition(unsigned n) const noexcept
{
    unsigned k {0};
    for (std::size_t i {0}; i < partitions_.size(); ++i) {
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

namespace debug {

template <typename S, typename T>
void print_alleles(S&& stream, const PartitionedGenotype<T>& genotype)
{
    for (std::size_t k {0}; k < genotype.num_partitions() - 1; ++k) {
        print_alleles(stream, genotype.partition(k));
        stream << " | ";
    }
    print_alleles(stream, genotype.partition(genotype.num_partitions() - 1));
}

template <typename S, typename T>
void print_variant_alleles(S&& stream, const PartitionedGenotype<T>& genotype)
{
    for (std::size_t k {0}; k < genotype.num_partitions() - 1; ++k) {
        print_variant_alleles(stream, genotype.partition(k));
        stream << " | ";
    }
    print_variant_alleles(stream, genotype.partition(genotype.num_partitions() - 1));
}

template <typename T>
void print_alleles(const PartitionedGenotype<T>& genotype)
{
    print_alleles(std::cout, genotype);
}

template <typename T>
void print_variant_alleles(const PartitionedGenotype<T>& genotype)
{
    print_variant_alleles(std::cout, genotype);
}

} // namespace debug

} // namespace octopus

#endif
