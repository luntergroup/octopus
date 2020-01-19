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
#include "haplotype.hpp"

namespace octopus {

template <typename IndexTp = std::size_t>
class IndexedHaplotype : public Mappable<IndexedHaplotype<IndexTp>>
                       , public Comparable<IndexedHaplotype<IndexTp>>
                       , public Indexed<IndexedHaplotype<IndexTp>>
{
public:
    using IndexType = IndexTp;
    
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

template <typename IndexType, typename H>
IndexedHaplotype<IndexType> index(H&& haplotype, IndexType index)
{
    return IndexedHaplotype<IndexType> {std::forward<H>(haplotype), index};
}

template <typename IndexType, typename InputIterator, typename OutputIterator,
          typename = std::enable_if_t<std::is_same<typename std::iterator_traits<InputIterator>::value_type, Haplotype>::value_type>>
OutputIterator
index(InputIterator first, InputIterator last, OutputIterator result, IndexType index = IndexType {0})
{
    return std::transform(first, last, result, [&] (auto&& h) { return index(std::forward<decltype(h)>(h), index++); });
}

template <typename IndexType = std::size_t>
std::vector<IndexedHaplotype<IndexType>>
index(const std::vector<Haplotype>& haplotypes, const IndexType index = IndexType {0})
{
    std::vector<IndexedHaplotype<IndexType>> result {};
    result.reserve(haplotypes.size());
    index(std::cbegin(haplotypes), std::cend(haplotypes),
          std::back_inserter(result), index);
    return result;
}

template <typename IndexType = std::size_t>
std::vector<IndexedHaplotype<IndexType>>
index(std::vector<Haplotype>&& haplotypes, const IndexType index = IndexType {0})
{
    std::vector<IndexedHaplotype<IndexType>> result {};
    result.reserve(haplotypes.size());
    index(std::make_move_iterator(std::begin(haplotypes)), std::make_move_iterator(std::end(haplotypes)),
          std::back_inserter(result), index);
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
unindex(const std::vector<Haplotype>& haplotypes)
{
    std::vector<Haplotype> result {};
    result.reserve(haplotypes.size());
    unindex(std::cbegin(haplotypes), std::cend(haplotypes),
            std::back_inserter(result));
    return result;
}

template <typename IndexType = std::size_t>
std::vector<IndexedHaplotype<IndexType>>
index(std::vector<Haplotype>&& haplotypes)
{
    std::vector<Haplotype> result {};
    result.reserve(haplotypes.size());
    unindex(std::make_move_iterator(std::begin(haplotypes)), std::make_move_iterator(std::end(haplotypes)),
            std::back_inserter(result));
    return result;
}

namespace detail {

template <typename>
struct IsIndexedHaplotype : std::false_type {};
template <typename IndexType>
struct IsIndexedHaplotype<IndexedHaplotype<IndexType>> : std::true_type {
    static_assert(is_indexed<IndexedHaplotype<IndexType>>(), "");
};

} // namespace detail

template <typename T>
constexpr bool is_indexed_haplotype = detail::IsIndexedHaplotype<T>::value;

} // namespace octopus

#endif
