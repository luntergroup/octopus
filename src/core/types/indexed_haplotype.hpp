// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef indexed_haplotype_hpp
#define indexed_haplotype_hpp

#include <cstddef>
#include <utility>
#include <type_traits>
#include <iterator>
#include <vector>
#include <functional>

#include "concepts/mappable.hpp"
#include "concepts/comparable.hpp"
#include "concepts/indexed.hpp"
#include "containers/mappable_block.hpp"
#include "haplotype.hpp"

namespace octopus {

template <typename IndexTp = unsigned,
          typename UseIndexEquality = std::true_type,
          typename IsSortedIndex = std::true_type>
class IndexedHaplotype : public Mappable<IndexedHaplotype<IndexTp>>
                       , public Comparable<IndexedHaplotype<IndexTp>>
                       , public Indexed<IndexedHaplotype<IndexTp>>
{
public:
    using IndexType = IndexTp;
    using IndexEquality = UseIndexEquality;
    using SortedIndex = IsSortedIndex;
    using MappingDomain = Haplotype::MappingDomain;
    
    IndexedHaplotype() = delete;
    
    IndexedHaplotype(const Haplotype& haplotype, IndexType index) noexcept
    : haplotype_ {haplotype}
    , index_ {index}
    {}
    
    IndexedHaplotype(const IndexedHaplotype&)            = default;
    IndexedHaplotype& operator=(const IndexedHaplotype&) = default;
    IndexedHaplotype(IndexedHaplotype&&)                 = default;
    IndexedHaplotype& operator=(IndexedHaplotype&&)      = default;
    
    ~IndexedHaplotype() = default;
    
    const Haplotype& haplotype() const noexcept { return haplotype_; }
    
    operator const Haplotype&() const noexcept { return haplotype(); }
    
    IndexType index() const noexcept { return index_; }
    
    decltype(auto) mapped_region() const noexcept { return haplotype().mapped_region(); }
    
    decltype(auto) contains(const ContigAllele& allele) const { return haplotype().contains(allele); }
    decltype(auto) contains(const Allele& allele) const  { return haplotype().contains(allele); }
    
    decltype(auto) includes(const ContigAllele& allele) const { return haplotype().includes(allele); }
    decltype(auto) includes(const Allele& allele) const  { return haplotype().includes(allele); }
    
    decltype(auto) sequence(const ContigRegion& region) const { return haplotype().sequence(region); }
    decltype(auto) sequence(const GenomicRegion& region) const { return haplotype().sequence(region); }
    decltype(auto) sequence() const noexcept { return haplotype().sequence(); }
    
    decltype(auto) sequence_size(const ContigRegion& region) const { return haplotype().sequence_size(region); }
    decltype(auto) sequence_size(const GenomicRegion& region) const { return haplotype().sequence_size(region); }
    
    decltype(auto) difference(const IndexedHaplotype& other) const { return haplotype().difference(other.haplotype()); }
    decltype(auto) cigar() const { return haplotype().cigar(); }
    
private:
    std::reference_wrapper<const Haplotype> haplotype_;
    IndexType index_;
};

template <typename IndexType>
bool operator==(const IndexedHaplotype<IndexType, std::true_type>& lhs,
                const IndexedHaplotype<IndexType, std::true_type>& rhs) noexcept
{
    return index_of(lhs) == index_of(rhs);
}
template <typename IndexType>
bool operator==(const IndexedHaplotype<IndexType, std::false_type>& lhs,
                const IndexedHaplotype<IndexType, std::false_type>& rhs) noexcept
{
    return lhs.haplotype() == rhs.haplotype();
}
template <typename IndexType, typename EqualityPolicy>
bool operator<(const IndexedHaplotype<IndexType, EqualityPolicy, std::true_type>& lhs,
               const IndexedHaplotype<IndexType, EqualityPolicy, std::true_type>& rhs) noexcept
{
    return index_of(lhs) < index_of(rhs);
}
template <typename IndexType, typename EqualityPolicy>
bool operator<(const IndexedHaplotype<IndexType, EqualityPolicy, std::false_type>& lhs,
               const IndexedHaplotype<IndexType, EqualityPolicy, std::false_type>& rhs)
{
    return lhs.haplotype() < rhs.haplotype();
}

template <typename IndexType>
IndexedHaplotype<IndexType> index(const Haplotype& haplotype, IndexType idx) noexcept
{
    return IndexedHaplotype<IndexType> {haplotype, idx};
}

template <typename IndexType, typename InputIterator, typename OutputIterator>
OutputIterator
index(InputIterator first, InputIterator last, OutputIterator result, IndexType init = IndexType {0}) noexcept
{
    return std::transform(first, last, result, [&] (const auto& h) { return index(h, init++); });
}

template <typename IndexType = IndexedHaplotype<>::IndexType>
std::vector<IndexedHaplotype<IndexType>>
index(const std::vector<Haplotype>& haplotypes, const IndexType init = IndexType {0})
{
    std::vector<IndexedHaplotype<IndexType>> result {};
    result.reserve(haplotypes.size());
    index(std::cbegin(haplotypes), std::cend(haplotypes), std::back_inserter(result), init);
    return result;
}

template <typename IndexType = IndexedHaplotype<>::IndexType>
MappableBlock<IndexedHaplotype<IndexType>>
index(const MappableBlock<Haplotype>& haplotypes, const IndexType init = IndexType {0})
{
    MappableBlock<IndexedHaplotype<IndexType>> result {mapped_region(haplotypes)};
    result.reserve(haplotypes.size());
    index(std::cbegin(haplotypes), std::cend(haplotypes), std::back_inserter(result), init);
    return result;
}

template <typename>
struct is_indexed_haplotype : std::false_type {};
template <typename IndexType>
struct is_indexed_haplotype<IndexedHaplotype<IndexType>> : std::true_type {
    static_assert(is_indexed_v<IndexedHaplotype<IndexType>>, "");
};

template <typename T>
constexpr bool is_indexed_haplotype_v = is_indexed_haplotype<T>::value;

template <typename IndexType, typename EqualityPolicy, typename IsSortedIndex, typename MappableType>
bool contains(const IndexedHaplotype<IndexType, EqualityPolicy, IsSortedIndex>& haplotype, const MappableType& other)
{
    return contains(haplotype.haplotype(), other);
}
template <typename IndexType, typename EqualityPolicy, typename IsSortedIndex, typename MappableType>
bool includes(const IndexedHaplotype<IndexType, EqualityPolicy, IsSortedIndex>& haplotype, const MappableType& other)
{
    return includes(haplotype.haplotype(), other);
}

} // namespace octopus

namespace std {

template <typename IndexType, typename EqualityPolicy, typename IsSortedIndex>
struct hash<octopus::IndexedHaplotype<IndexType, EqualityPolicy, IsSortedIndex>>
{
    size_t operator()(const octopus::IndexedHaplotype<IndexType, EqualityPolicy, IsSortedIndex>& haplotpe) const noexcept
    {
        return haplotpe.index();
    }
};

} // namespace std

namespace boost {

template <typename IndexType, typename EqualityPolicy, typename IsSortedIndex>
struct hash<octopus::IndexedHaplotype<IndexType, EqualityPolicy, IsSortedIndex>>
{
    std::size_t operator()(const octopus::IndexedHaplotype<IndexType, EqualityPolicy, IsSortedIndex>& h) const noexcept
    {
        return std::hash<octopus::IndexedHaplotype<IndexType, EqualityPolicy, IsSortedIndex>>()(h);
    }
};

} // namespace boost

#endif
