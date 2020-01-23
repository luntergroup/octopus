// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef indexed_haplotype_hpp
#define indexed_haplotype_hpp

#include <cstddef>
#include <utility>
#include <type_traits>
#include <iterator>
#include <vector>

#include "concepts/mappable.hpp"
#include "concepts/comparable.hpp"
#include "concepts/indexed.hpp"
#include "containers/mappable_block.hpp"
#include "haplotype.hpp"

namespace octopus {

template <typename IndexTp = std::size_t,
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
    
    template <typename H>
    IndexedHaplotype(H&& haplotype, IndexType index)
    : haplotype_ {std::forward<H>(haplotype)}
    , index_ {index}
    {}
    
    IndexedHaplotype(const IndexedHaplotype&)            = default;
    IndexedHaplotype& operator=(const IndexedHaplotype&) = default;
    IndexedHaplotype(IndexedHaplotype&&)                 = default;
    IndexedHaplotype& operator=(IndexedHaplotype&&)      = default;
    
    ~IndexedHaplotype() = default;
    
    decltype(auto) mapped_region() const noexcept { return haplotype_.mapped_region(); }
    
    IndexType index() const noexcept { return index_; }
    
    Haplotype& haplotype() noexcept { return haplotype_; }
    const Haplotype& haplotype() const noexcept { return haplotype_; }
    
    operator Haplotype&() noexcept { return haplotype_; }
    operator const Haplotype&() const noexcept { return haplotype_; }
    
    decltype(auto) contains(const ContigAllele& allele) const { return haplotype_.contains(allele); }
    decltype(auto) contains(const Allele& allele) const  { return haplotype_.contains(allele); }
    
    decltype(auto) includes(const ContigAllele& allele) const { return haplotype_.includes(allele); }
    decltype(auto) includes(const Allele& allele) const  { return haplotype_.includes(allele); }
    
    decltype(auto) sequence(const ContigRegion& region) const { return haplotype_.sequence(region); }
    decltype(auto) sequence(const GenomicRegion& region) const { return haplotype_.sequence(region); }
    decltype(auto) sequence() const noexcept { return haplotype_.sequence(); }
    
    decltype(auto) sequence_size(const ContigRegion& region) const { return haplotype_.sequence_size(region); }
    decltype(auto) sequence_size(const GenomicRegion& region) const { return haplotype_.sequence_size(region); }
    
    decltype(auto) difference(const IndexedHaplotype& other) const { return haplotype_.difference(other); }
    decltype(auto) cigar() const { return haplotype_.cigar(); }
    
private:
    Haplotype haplotype_;
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
                const IndexedHaplotype<IndexType, std::false_type>& rhs)
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

template <typename IndexType, typename H>
IndexedHaplotype<IndexType> index(H&& haplotype, IndexType idx)
{
    return IndexedHaplotype<IndexType> {std::forward<H>(haplotype), idx};
}

template <typename IndexType, typename InputIterator, typename OutputIterator>
OutputIterator
index(InputIterator first, InputIterator last, OutputIterator result, IndexType init = IndexType {0})
{
    return std::transform(first, last, result, [&] (auto&& h) { return index(std::forward<decltype(h)>(h), init++); });
}

template <typename IndexType = std::size_t>
std::vector<IndexedHaplotype<IndexType>>
index(const std::vector<Haplotype>& haplotypes, const IndexType init = IndexType {0})
{
    std::vector<IndexedHaplotype<IndexType>> result {};
    result.reserve(haplotypes.size());
    index(std::cbegin(haplotypes), std::cend(haplotypes),
          std::back_inserter(result), init);
    return result;
}

template <typename IndexType = std::size_t>
std::vector<IndexedHaplotype<IndexType>>
index(std::vector<Haplotype>&& haplotypes, const IndexType init = IndexType {0})
{
    std::vector<IndexedHaplotype<IndexType>> result {};
    result.reserve(haplotypes.size());
    index(std::make_move_iterator(std::begin(haplotypes)), std::make_move_iterator(std::end(haplotypes)),
          std::back_inserter(result), init);
    return result;
}

template <typename IndexType = std::size_t>
MappableBlock<IndexedHaplotype<IndexType>>
index(const MappableBlock<Haplotype>& haplotypes, const IndexType init = IndexType {0})
{
    MappableBlock<IndexedHaplotype<IndexType>> result {mapped_region(haplotypes)};
    result.reserve(haplotypes.size());
    index(std::cbegin(haplotypes), std::cend(haplotypes),
          std::back_inserter(result), init);
    return result;
}

template <typename IndexType = std::size_t>
MappableBlock<IndexedHaplotype<IndexType>>
index(MappableBlock<Haplotype>&& haplotypes, const IndexType init = IndexType {0})
{
    MappableBlock<IndexedHaplotype<IndexType>> result {mapped_region(haplotypes)};
    result.reserve(haplotypes.size());
    index(std::make_move_iterator(std::begin(haplotypes)), std::make_move_iterator(std::end(haplotypes)),
          std::back_inserter(result), init);
    return result;
}

template <typename IndexType>
auto unindex(const IndexedHaplotype<IndexType>& haplotype)
{
    return haplotype.haplotype();
}
template <typename IndexType>
auto unindex(IndexedHaplotype<IndexType>&& haplotype)
{
    return std::move(haplotype.haplotype());
}

template <typename InputIterator, typename OutputIterator>
OutputIterator
unindex(InputIterator first, InputIterator last, OutputIterator result)
{
    return std::transform(first, last, result, [&] (auto&& h) { return unindex(std::forward<decltype(h)>(h)); });
}

template <typename IndexType>
std::vector<Haplotype>
unindex(const std::vector<IndexedHaplotype<IndexType>>& haplotypes)
{
    std::vector<Haplotype> result {};
    result.reserve(haplotypes.size());
    unindex(std::cbegin(haplotypes), std::cend(haplotypes),
            std::back_inserter(result));
    return result;
}

template <typename IndexType = std::size_t>
std::vector<Haplotype>
unindex(std::vector<IndexedHaplotype<IndexType>>&& haplotypes)
{
    std::vector<Haplotype> result {};
    result.reserve(haplotypes.size());
    unindex(std::make_move_iterator(std::begin(haplotypes)), std::make_move_iterator(std::end(haplotypes)),
            std::back_inserter(result));
    return result;
}

template <typename IndexType>
MappableBlock<Haplotype>
unindex(const MappableBlock<IndexedHaplotype<IndexType>>& haplotypes)
{
    MappableBlock<Haplotype> result {};
    result.reserve(haplotypes.size());
    unindex(std::cbegin(haplotypes), std::cend(haplotypes),
            std::back_inserter(result));
    return result;
}

template <typename IndexType = std::size_t>
MappableBlock<Haplotype>
unindex(MappableBlock<IndexedHaplotype<IndexType>>&& haplotypes)
{
    MappableBlock<Haplotype> result {};
    result.reserve(haplotypes.size());
    unindex(std::make_move_iterator(std::begin(haplotypes)), std::make_move_iterator(std::end(haplotypes)),
            std::back_inserter(result));
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

} // namespace octopus

#endif
