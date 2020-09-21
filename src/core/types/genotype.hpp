// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef genotype_hpp
#define genotype_hpp

#include <vector>
#include <deque>
#include <memory>
#include <initializer_list>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <cstddef>
#include <iterator>
#include <algorithm>
#include <type_traits>
#include <ostream>
#include <cassert>
#include <iostream>

#include <boost/functional/hash.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "concepts/equitable.hpp"
#include "concepts/mappable.hpp"
#include "containers/mappable_block.hpp"
#include "utils/reorder.hpp"
#include "allele.hpp"
#include "haplotype.hpp"
#include "indexed_haplotype.hpp"
#include "shared_haplotype.hpp"

namespace octopus {

namespace detail {

template <typename T>
constexpr bool is_allele_v = std::is_same<T, Allele>::value || std::is_same<T, ContigAllele>::value;

template <typename MappaleType, typename = void>
struct is_haplotype_like : public std::false_type {};
template <>
struct is_haplotype_like<Haplotype> : std::true_type {};
template <>
struct is_haplotype_like<SharedHaplotype> : std::true_type {};
template <typename MappaleType>
struct is_haplotype_like<MappaleType, std::enable_if_t<is_indexed_haplotype_v<MappaleType>>> : std::true_type {};

template <typename T>
constexpr bool is_haplotype_like_v = is_haplotype_like<T>::value;

} // namespace detail

template <typename T>
constexpr bool is_genotypeable_v = detail::is_allele_v<T> || detail::is_haplotype_like_v<T>;

template <typename MappableType> class Genotype;

namespace detail {

template <typename T> Genotype<T> collapse(const Genotype<T>& genotype, std::true_type);

} // namespace detail

template <typename MappableType>
class Genotype : public Equitable<Genotype<MappableType>>, public Mappable<Genotype<MappableType>>
{
public:
    static_assert(is_genotypeable_v<MappableType>, "");
    
    using ElementType   = MappableType;
    using MappingDomain = RegionType<ElementType>;
    
    using ordered = detail::is_haplotype_like<MappableType>;
    
    using value_type = ElementType; // for use with genetic algorithms
    
    Genotype() = default;
    
    Genotype(unsigned ploidy);
    Genotype(unsigned ploidy, const MappableType& init);
    Genotype(std::initializer_list<MappableType> elements);
    template <typename InputIterator>
    Genotype(InputIterator first, InputIterator last);
    
    Genotype(const Genotype&)            = default;
    Genotype& operator=(const Genotype&) = default;
    Genotype(Genotype&&)                 = default;
    Genotype& operator=(Genotype&&)      = default;
    
    ~Genotype() = default;
    
    template <typename T> void emplace(T&& element);
    
    ElementType& operator[](unsigned n) noexcept;
    const ElementType& operator[](unsigned n) const noexcept;
    
    const MappingDomain& mapped_region() const noexcept;
    
    unsigned ploidy() const noexcept;
    
    void reorder(const std::vector<unsigned>& order) { reorder(order, ordered {}); }
    void collapse() { collapse(ordered {}); }
    
    template <typename T> friend Genotype<T> detail::collapse(const Genotype<T>& genotype, std::true_type);

    auto capacity() const noexcept { return elements_.capacity(); }
    void shrink_to_fit() { elements_.shrink_to_fit(); }
    
private:
    std::vector<MappableType> elements_;
    
    void init() { init(ordered{}); }
    void init(std::true_type);
    void init(std::false_type) const noexcept {};
    
    template <typename T> void emplace(T&& element, std::false_type);
    template <typename T> void emplace(T&& element, std::true_type);
    void reorder(const std::vector<unsigned>& order, std::false_type);
    void collapse(std::true_type);

public:
    using iterator       = typename decltype(elements_)::iterator;
    using const_iterator = typename decltype(elements_)::const_iterator;
    
    iterator begin() noexcept { return std::begin(elements_); }
    iterator end() noexcept { return std::end(elements_); }
    const_iterator begin() const noexcept { return std::cbegin(elements_); }
    const_iterator end() const noexcept { return std::cend(elements_); }
    const_iterator cbegin() const noexcept { return std::cbegin(elements_); }
    const_iterator cend() const noexcept { return std::cend(elements_); }
};

// Genotype<MappableType>

template <typename MappableType>
Genotype<MappableType>::Genotype(const unsigned ploidy)
: elements_ {}
{
    elements_.reserve(ploidy);
}

template <typename MappableType>
Genotype<MappableType>::Genotype(const unsigned ploidy, const MappableType& init)
: elements_ {ploidy, init}
{}

template <typename MappableType>
Genotype<MappableType>::Genotype(std::initializer_list<MappableType> elements)
: elements_ {elements}
{
    init();
}

template <typename MappableType>
template <typename InputIterator>
Genotype<MappableType>::Genotype(InputIterator first, InputIterator last)
: elements_ {first, last}
{
    init();
}

template <typename MappableType>
void Genotype<MappableType>::init(std::true_type)
{
    std::sort(std::begin(elements_), std::end(elements_));
}

template <typename MappableType>
const typename Genotype<MappableType>::MappingDomain& Genotype<MappableType>::mapped_region() const noexcept
{
    return octopus::mapped_region(elements_[0]);
}

template <typename MappableType>
MappableType& Genotype<MappableType>::operator[](const unsigned n) noexcept
{
    return elements_[n];
}
template <typename MappableType>
const MappableType& Genotype<MappableType>::operator[](const unsigned n) const noexcept
{
    return elements_[n];
}

template <typename MappableType>
unsigned Genotype<MappableType>::ploidy() const noexcept
{
    return elements_.size();
}

namespace detail {

template <typename T, typename Container, typename BinaryPredicate>
void insert_sorted(T&& value, Container& values, BinaryPredicate compare)
{
    auto position = std::upper_bound(std::cbegin(values), std::cend(values), value, compare);
    values.insert(position, std::forward<T>(value));
}

} // namespace detail

template <typename MappableType>
template <typename T>
void Genotype<MappableType>::emplace(T&& element, std::false_type)
{
    elements_.emplace_back(std::forward<T>(element));
}
template <typename MappableType>
template <typename T>
void Genotype<MappableType>::emplace(T&& element, std::true_type)
{
    detail::insert_sorted(std::forward<T>(element), elements_, std::less<> {});
}
template <typename MappableType>
template <typename T>
void Genotype<MappableType>::emplace(T&& element)
{
    emplace(std::forward<T>(element), ordered {});
}

template <typename MappableType>
void Genotype<MappableType>::reorder(const std::vector<unsigned>& order, std::false_type)
{
    assert(order.size() == elements_.size());
    octopus::reorder(std::cbegin(order), std::cend(order), std::begin(elements_));
}

template <typename MappableType>
void Genotype<MappableType>::collapse(std::true_type)
{
    elements_.erase(std::unique(std::begin(elements_), std::end(elements_)), std::end(elements_));
}

// free functions

template <typename MappableType>
auto ploidy(const Genotype<MappableType>& genotype) noexcept
{
    return genotype.ploidy();
}

template <typename MappableType>
bool is_max_zygosity(const Genotype<MappableType>& genotype)
{
    return zygosity(genotype) == ploidy(genotype);
}

namespace detail {

template <typename MappableType>
bool is_homozygous(const Genotype<MappableType>& genotype, std::false_type)
{
    return std::adjacent_find(std::cbegin(genotype), std::cend(genotype), std::not_equal_to<MappableType>{}) == std::cend(genotype);
}
template <typename MappableType>
bool is_homozygous(const Genotype<MappableType>& genotype, std::true_type)
{
    return genotype[0] == genotype[genotype.ploidy() - 1];
}

template <typename MappableType>
unsigned zygosity(const Genotype<MappableType>& genotype, std::false_type)
{
    if (genotype.ploidy() == 1 || is_homozygous(genotype, std::false_type {})) {
        return 1;
    } else if (genotype.ploidy() == 2) {
        return 2;
    } else if (genotype.ploidy == 3) {
        if (genotype[0] == genotype[1] || genotype[0] == genotype[2]) {
            return 2;
        } else {
            return 3;
        }
    }
    return copy_unique(genotype, std::false_type {}).size();
}
template <typename MappableType>
unsigned zygosity(const Genotype<MappableType>& genotype, std::true_type)
{
    if (genotype.ploidy() < 2) return genotype.ploidy();
    unsigned result {0};
    for (auto element_itr = std::cbegin(genotype); element_itr != std::cend(genotype); ++result) {
        // naive algorithm faster in practice than binary searching
        const auto is_same_element = [=] (const MappableType& element) { return element == *element_itr; };
        element_itr = std::find_if_not(std::next(element_itr), std::cend(genotype), is_same_element);
    }
    return result;
}

template <typename MappableType>
bool contains(const Genotype<MappableType>& genotype, const MappableType& element, std::false_type)
{
    return std::find(std::cbegin(genotype), std::cend(genotype), element) != std::cend(genotype);
}
template <typename MappableType>
bool contains(const Genotype<MappableType>& genotype, const MappableType& element, std::true_type)
{
    return std::binary_search(std::cbegin(genotype), std::cend(genotype), element);
}

template <typename MappableType>
unsigned count(const Genotype<MappableType>& genotype, const MappableType& element, std::false_type)
{
    return std::count(std::cbegin(genotype), std::cend(genotype), element);
}
template <typename MappableType>
unsigned count(const Genotype<MappableType>& genotype, const MappableType& element, std::true_type)
{
    const auto equal_range = std::equal_range(std::cbegin(genotype), std::cend(genotype), element);
    return std::distance(equal_range.first, equal_range.second);
}

//template <typename MappableType>
//Genotype<MappableType> collapse(const Genotype<MappableType>& genotype, std::false_type)
//{
//    return {};
//}
template <typename MappableType>
Genotype<MappableType> collapse(const Genotype<MappableType>& genotype, std::true_type)
{
    if (is_max_zygosity(genotype)) return genotype;
    Genotype<MappableType> result {};
    result.elements_.reserve(ploidy(genotype));
    std::unique_copy(std::cbegin(genotype), std::cend(genotype), std::back_inserter(result.elements_));
    result.elements_.shrink_to_fit();
    return result;
}

//template <typename MappableType>
//std::vector<unsigned> unique_counts(const Genotype<MappableType>& genotype, std::false_type)
//{
//    std::vector<unsigned> result {};
//    return result;
//}
template <typename MappableType>
std::vector<unsigned> unique_counts(const Genotype<MappableType>& genotype, std::true_type)
{
    std::vector<unsigned> result {};
    result.reserve(genotype.ploidy());
    for (auto element_itr = std::cbegin(genotype); element_itr != std::cend(genotype);) {
        const auto is_same_element = [=] (const auto& element) { return element == *element_itr; };
        const auto next_element_itr = std::find_if_not(std::next(element_itr), std::cend(genotype), is_same_element);
        result.push_back(std::distance(element_itr, next_element_itr));
        element_itr = next_element_itr;
    }
    return result;
}

template <typename InputIt1, typename InputIt2,
          typename Compare>
bool set_equal(InputIt1 first1, InputIt1 last1, 
               InputIt2 first2, InputIt2 last2,
               Compare cmp)
{
    while (first1 != last1) {
        if (first2 == last2) return false;
        if (!cmp(*first1, *first2)) return false;
        const auto& match = *first1;
        const auto is_duplicate = [&] (const auto& v) { return cmp(v, match); };
        first1 = std::find_if_not(std::next(first1), last1, is_duplicate);
        first2 = std::find_if_not(std::next(first2), last2, is_duplicate);
    }
    return first2 == last2;
}

template <typename InputIt1, typename InputIt2>
bool set_equal(InputIt1 first1, InputIt1 last1, 
               InputIt2 first2, InputIt2 last2)
{
    return set_equal(first1, last1, first2, last2, std::equal_to<> {});
}

template <typename MappableType>
bool have_same_elements(const Genotype<MappableType>& lhs, const Genotype<MappableType>& rhs, std::true_type)
{
    return set_equal(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs), std::cend(rhs));
}

} // namespace detail

template <typename MappableType>
using is_ordered = typename Genotype<MappableType>::ordered;

template <typename MappableType>
constexpr bool is_ordered_v = is_ordered<MappableType>::value;

template <typename MappableType>
bool is_homozygous(const Genotype<MappableType>& genotype)
{
    return detail::is_homozygous(genotype, is_ordered<MappableType> {});
}

template <typename MappableType>
bool is_heterozygous(const Genotype<MappableType>& genotype)
{
    return !is_homozygous(genotype);
}

template <typename T, std::enable_if_t<detail::is_haplotype_like_v<T>, int> = 0>
bool is_homozygous_reference(const Genotype<T>& genotype)
{
    return is_homozygous(genotype) && is_reference(genotype[0]);
}

template <typename MappableType>
unsigned zygosity(const Genotype<MappableType>& genotype)
{
    return detail::zygosity(genotype, is_ordered<MappableType> {});
}

template <typename MappableType>
bool contains(const Genotype<MappableType>& genotype, const MappableType& element)
{
    return detail::contains(genotype, element, is_ordered<MappableType> {});
}

template <typename MappableType>
unsigned count(const Genotype<MappableType>& genotype, const MappableType& element)
{
    return detail::count(genotype, element, is_ordered<MappableType> {});
}

template <typename MappableType>
Genotype<MappableType> collapse(const Genotype<MappableType>& genotype)
{
    return detail::collapse(genotype, is_ordered<MappableType> {});
}

template <typename MappableType>
void collapse(Genotype<MappableType>& genotype)
{
    genotype.collapse();
}

template <typename MappableType>
std::vector<unsigned> unique_counts(const Genotype<MappableType>& genotype)
{
    return detail::unique_counts(genotype, is_ordered<MappableType> {});
}

template <typename MappableType>
bool have_same_elements(const Genotype<MappableType>& lhs, const Genotype<MappableType>& rhs)
{
    return detail::have_same_elements(lhs, rhs, is_ordered<MappableType> {});
}

template <typename MappableType>
bool is_haploid(const Genotype<MappableType>& genotype) noexcept
{
    return genotype.ploidy() == 1;
}
template <typename MappableType>
bool is_diploid(const Genotype<MappableType>& genotype) noexcept
{
    return genotype.ploidy() == 2;
}
template <typename MappableType>
bool is_triploid(const Genotype<MappableType>& genotype) noexcept
{
    return genotype.ploidy() == 3;
}
template <typename MappableType>
bool is_tetraploid(const Genotype<MappableType>& genotype) noexcept
{
    return genotype.ploidy() == 4;
}

namespace detail {

template <typename MappableType, typename UnaryFunction>
auto transform(const Genotype<MappableType>& genotype, UnaryFunction&& f, std::false_type)
{
    Genotype<std::result_of_t<UnaryFunction(const MappableType&)>> result {genotype.ploidy()};
    for (const auto& element : genotype) {
        result.emplace(f(element));
    }
    return result;
}
template <typename MappableType, typename UnaryFunction>
auto transform(const Genotype<MappableType>& genotype, UnaryFunction&& f, std::true_type)
{
    Genotype<std::result_of_t<UnaryFunction(const MappableType&)>> result {genotype.ploidy()};
    for (auto element_itr = std::cbegin(genotype); element_itr != std::cend(genotype);) {
        auto transformed_element = f(*element_itr);
        ++element_itr;
        for (; element_itr != std::cend(genotype) && *element_itr == *std::prev(element_itr); ++element_itr) {
            result.emplace(transformed_element);
        }
        result.emplace(std::move(transformed_element));
    }
    return result;
}

} // namespace detail

template <typename MappableType, typename UnaryFunction>
auto transform(const Genotype<MappableType>& genotype, UnaryFunction&& f)
{
    return detail::transform(genotype, std::forward<UnaryFunction>(f), is_ordered<MappableType> {});
}

template <typename MappableType>
const Genotype<MappableType>& genotype_cast(const Genotype<MappableType>& genotype) noexcept
{
    return genotype;
}
template <typename MappableType2, typename MappableType1>
Genotype<MappableType2> genotype_cast(const Genotype<MappableType1>& genotype)
{
    return transform(genotype, [] (const MappableType1& element) { return static_cast<MappableType2>(element); });
}

template <typename MappableType, typename MappableType_>
Genotype<MappableType> copy(const Genotype<MappableType_>& genotype, const GenomicRegion& region)
{
    static_assert(is_ordered_v<MappableType_>, "");
    return transform(genotype, [&] (const auto& element) { return copy<MappableType>(element, region); });
}

namespace detail {

template <typename ForwardIterator,
          typename UnaryFunction1,
          typename UnaryFunction2,
          typename OutputIterator>
auto transform_each_naive(ForwardIterator first, ForwardIterator last,
                          UnaryFunction1&& element_f,
                          UnaryFunction2&& genotype_f,
                          OutputIterator result)
{
    return std::transform(first, last, result, [&] (const auto& genotype) { return genotype_f(transform(genotype, element_f)); });
}

template <typename MappableType, typename UnaryFunction>
using MapCache = std::unordered_map<MappableType, std::result_of_t<UnaryFunction(const MappableType&)>>;
template <typename MappableType, typename UnaryFunction>
using LinearCache = std::vector<boost::optional<std::result_of_t<UnaryFunction(const MappableType&)>>>;

template <typename ElementType, typename UnaryFunction>
auto init_cache(std::size_t size_hint, std::true_type) { return LinearCache<ElementType, UnaryFunction>(size_hint); }
template <typename ElementType, typename UnaryFunction>
auto init_cache(std::size_t size_hint, std::false_type) { return MapCache<ElementType, UnaryFunction>(size_hint); }
template <typename ElementType, typename UnaryFunction>
auto init_cache(std::size_t size_hint = 100)
{
    return init_cache<ElementType, UnaryFunction>(size_hint, is_indexed<ElementType> {});
}

template <typename MappableType, typename UnaryFunction>
decltype(auto)
copy_from_cache(const MappableType& element, UnaryFunction&& f, MapCache<MappableType, UnaryFunction>& cache)
{
    const auto itr = cache.find(element);
    if (itr == std::cend(cache)) {
        auto result = f(element);
        cache.emplace(std::piecewise_construct, std::forward_as_tuple(element), std::forward_as_tuple(result));
        return result;
    } else {
        return itr->second;
    }
}
template <typename IndexedType, typename UnaryFunction>
decltype(auto)
copy_from_cache(const IndexedType& element, UnaryFunction&& f, LinearCache<IndexedType, UnaryFunction>& cache)
{
    if (index_of(element) >= cache.size()) {
        cache.resize(2 * (index_of(element) + 1));
    }
    if (!cache[index_of(element)]) {
        cache[index_of(element)] = f(element);
    }
    return *cache[index_of(element)];
}

template <typename MappableType, typename UnaryFunction, typename CacheType>
auto
transform(const Genotype<MappableType>& genotype, UnaryFunction&& f, CacheType& cache)
{
    static_assert(is_ordered_v<MappableType>, "");
    Genotype<std::result_of_t<UnaryFunction(const MappableType&)>> result {genotype.ploidy()};
    for (auto element_itr = std::cbegin(genotype); element_itr != std::cend(genotype);) {
        auto element_copy = copy_from_cache(*element_itr, f, cache);
        ++element_itr;
        for (; element_itr != std::cend(genotype) && *element_itr == *std::prev(element_itr); ++element_itr) {
            result.emplace(element_copy);
        }
        result.emplace(std::move(element_copy));
    }
    return result;
}

template <typename ForwardIterator,
          typename UnaryFunction1,
          typename UnaryFunction2,
          typename OutputIterator>
auto transform_each_cached(ForwardIterator first, ForwardIterator last,
                          UnaryFunction1&& element_f,
                          UnaryFunction2&& genotype_f,
                          OutputIterator result)
{
    if (first == last) return result;
    using ElementType = typename std::iterator_traits<ForwardIterator>::value_type::value_type;
    const auto cache_hint = std::min(static_cast<std::size_t>(std::distance(first, last)), std::size_t {100});
    auto cache = init_cache<ElementType, UnaryFunction1>(cache_hint);
    return std::transform(first, last, result, [&] (const auto& genotype) { return genotype_f(transform(genotype, element_f, cache)); });
}

} // namespace detail

template <typename ForwardIterator,
          typename UnaryFunction1,
          typename UnaryFunction2,
          typename OutputIterator>
auto transform_each(ForwardIterator first, ForwardIterator last,
                    UnaryFunction1&& element_f, 
                    UnaryFunction2&& genotype_f,
                    OutputIterator result)
{
    if (std::distance(first, last) < 5) {
        return detail::transform_each_naive(first, last, std::forward<UnaryFunction1>(element_f), std::forward<UnaryFunction2>(genotype_f), result);
    } else {
        return detail::transform_each_cached(first, last, std::forward<UnaryFunction1>(element_f), std::forward<UnaryFunction2>(genotype_f), result);
    }
}

template <typename ForwardIterator, typename UnaryFunction, typename OutputIterator>
auto transform_each(ForwardIterator first, ForwardIterator last, UnaryFunction&& element_f, OutputIterator result)
{
    return transform_each(first, last, std::forward<UnaryFunction>(element_f), [] (auto g) { return g; }, result);
}

template <typename Range, typename UnaryFunction, typename OutputIterator>
auto transform_each(const Range& genotypes, UnaryFunction&& f, OutputIterator result)
{
    return transform_each(std::cbegin(genotypes), std::cend(genotypes), std::forward<UnaryFunction>(f), result);
}

template <typename MappableType, typename MappableType_>
auto copy_each(const std::vector<Genotype<MappableType_>>& genotypes, const GenomicRegion& region)
{
    std::vector<Genotype<MappableType>> result {};
    result.reserve(genotypes.size());
    transform_each(genotypes, [&] (const auto& element) { return copy<MappableType>(element, region); }, std::back_inserter(result));
    return result;
}
template <typename MappableType, typename MappableType_>
auto copy_each(const MappableBlock<Genotype<MappableType_>>& genotypes, const GenomicRegion& region)
{
    MappableBlock<Genotype<MappableType>> result {region};
    result.reserve(genotypes.size());
    transform_each(genotypes, [&] (const auto& element) { return copy<MappableType>(element, region); }, std::back_inserter(result));
    return result;
}

template <typename MappableType, typename MappableType_>
Genotype<MappableType> copy(const Genotype<MappableType_>& genotype, const std::vector<GenomicRegion>& regions)
{
    static_assert(is_ordered_v<MappableType_>, "");
    return transform(genotype, [&] (const auto& element) { return copy<MappableType>(element, regions); });
}

template <typename MappableType, typename MappableType_>
auto copy_each(const std::vector<Genotype<MappableType_>>& genotypes, const std::vector<GenomicRegion>& regions)
{
    std::vector<Genotype<MappableType>> result {};
    result.reserve(genotypes.size());
    transform_each(genotypes, [&] (const auto& element) -> MappableType { return copy(element, regions); }, std::back_inserter(result));
    return result;
}
template <typename MappableType, typename MappableType_>
auto copy_each(const MappableBlock<Genotype<MappableType_>>& genotypes, const std::vector<GenomicRegion>& regions)
{
    auto copy_region = genotypes.empty() ? mapped_region(genotypes) : closed_region(regions.front(), regions.back());
    MappableBlock<Genotype<MappableType>> result {std::move(copy_region)};
    result.reserve(genotypes.size());
    transform_each(genotypes, [&] (const auto& element) -> MappableType { return copy(element, regions); }, std::back_inserter(result));
    return result;
}

namespace detail {

template <typename MappableType, typename UnaryPredicate>
bool any_of(const Genotype<MappableType>& genotype, UnaryPredicate&& pred, std::true_type)
{
    for (auto element_itr = std::cbegin(genotype); element_itr != std::cend(genotype);) {
        if (pred(*element_itr)) return true;
        const auto is_same_element = [element_itr] (const auto& x) { return x == *element_itr; };
        element_itr = std::find_if_not(std::next(element_itr), std::cend(genotype), is_same_element);
    }
    return false;
}
template <typename MappableType, typename UnaryPredicate>
bool any_of(const Genotype<MappableType>& genotype, UnaryPredicate&& pred, std::false_type)
{
    return std::any_of(std::cbegin(genotype), std::cend(genotype), std::forward<UnaryPredicate>(pred));
}

template <typename MappableType, typename UnaryPredicate>
bool all_of(const Genotype<MappableType>& genotype, UnaryPredicate&& pred, std::true_type)
{
    for (auto element_itr = std::cbegin(genotype); element_itr != std::cend(genotype);) {
        if (!pred(*element_itr)) return false;
        const auto is_same_element = [element_itr] (const auto& x) { return x == *element_itr; };
        element_itr = std::find_if_not(std::next(element_itr), std::cend(genotype), is_same_element);
    }
    return true;
}
template <typename MappableType, typename UnaryPredicate>
bool all_of(const Genotype<MappableType>& genotype, UnaryPredicate&& pred, std::false_type)
{
    return std::all_of(std::cbegin(genotype), std::cend(genotype), std::forward<UnaryPredicate>(pred));
}

template <typename MappableType, typename UnaryPredicate, typename CacheType>
unsigned count_if(const Genotype<MappableType>& genotype, UnaryPredicate&& pred, std::true_type)
{
    unsigned result {0};
    for (auto element_itr = std::cbegin(genotype); element_itr != std::cend(genotype);) {
        if (pred(*element_itr++)) {
            ++result;
            for (; element_itr != std::cend(genotype) && *element_itr == *std::prev(element_itr); ++element_itr, ++result);
        } else {
            for (; element_itr != std::cend(genotype) && *element_itr == *std::prev(element_itr); ++element_itr);
        }
    }
    return result;
}

template <typename MappableType, typename UnaryPredicate, typename CacheType>
unsigned count_if(const Genotype<MappableType>& genotype, UnaryPredicate&& pred, std::false_type)
{
    return std::count_if(std::cbegin(genotype), std::cend(genotype), std::forward<UnaryPredicate>(pred));
}

} // namespace detail

template <typename MappableType, typename UnaryPredicate>
bool any_of(const Genotype<MappableType>& genotype, UnaryPredicate&& pred)
{
    return detail::any_of(genotype, std::forward<UnaryPredicate>(pred), is_ordered<MappableType> {});
}
template <typename MappableType, typename UnaryPredicate>
bool all_of(const Genotype<MappableType>& genotype, UnaryPredicate&& pred)
{
    return detail::all_of(genotype, std::forward<UnaryPredicate>(pred), is_ordered<MappableType> {});
}

template <typename T, std::enable_if_t<detail::is_haplotype_like_v<T>, int> = 0>
bool contains(const Genotype<T>& genotype, const Allele& allele)
{
    return any_of(genotype, [&] (const T& haplotype) { return haplotype.contains(allele); });
}

template <typename T, std::enable_if_t<detail::is_haplotype_like_v<T>, int> = 0>
bool includes(const Genotype<T>& genotype, const Allele& allele)
{
    return any_of(genotype, [&] (const T& haplotype) { return haplotype.includes(allele); });
}

namespace detail {

template <typename MappableType, typename UnaryPredicate, typename CacheType>
bool any_of(const Genotype<MappableType>& genotype, UnaryPredicate&& pred, CacheType& cache)
{
    static_assert(is_ordered_v<MappableType>, "");
    for (auto element_itr = std::cbegin(genotype); element_itr != std::cend(genotype);) {
        if (copy_from_cache(*element_itr, pred, cache)) return true;
        ++element_itr;
        for (; element_itr != std::cend(genotype) && *element_itr == *std::prev(element_itr); ++element_itr);
    }
    return false;
}

template <typename MappableType, typename UnaryPredicate, typename CacheType>
bool all_of(const Genotype<MappableType>& genotype, UnaryPredicate&& pred, CacheType& cache)
{
    static_assert(is_ordered_v<MappableType>, "");
    for (auto element_itr = std::cbegin(genotype); element_itr != std::cend(genotype);) {
        if (!copy_from_cache(*element_itr, pred, cache)) return false;
        ++element_itr;
        for (; element_itr != std::cend(genotype) && *element_itr == *std::prev(element_itr); ++element_itr);
    }
    return true;
}

template <typename MappableType, typename UnaryPredicate, typename CacheType>
unsigned count_if(const Genotype<MappableType>& genotype, UnaryPredicate&& pred, CacheType& cache)
{
    static_assert(is_ordered_v<MappableType>, "");
    unsigned result {0};
    for (auto element_itr = std::cbegin(genotype); element_itr != std::cend(genotype);) {
        if (copy_from_cache(*element_itr++, pred, cache)) {
            ++result;
            for (; element_itr != std::cend(genotype) && *element_itr == *std::prev(element_itr); ++element_itr, ++result);
        } else {
            for (; element_itr != std::cend(genotype) && *element_itr == *std::prev(element_itr); ++element_itr);
        }
    }
    return result;
}

} // namespace detail

template <typename InputIterator, typename UnaryPredicate, typename BinaryFunction, typename UnaryFunction>
void for_each_any_of(InputIterator first, InputIterator last,
                     UnaryPredicate&& pred,
                     BinaryFunction&& visitor,
                     UnaryFunction&& adaptor)
{
    using V = typename std::iterator_traits<InputIterator>::value_type;
    using InputGenotypeType = std::result_of_t<UnaryFunction(V)>;
    using InputElementType = typename InputGenotypeType::ElementType;
    static_assert(detail::is_haplotype_like_v<InputElementType>, "");
    if (first == last) return;
    const auto cache_hint = std::min(static_cast<std::size_t>(std::distance(first, last)), std::size_t {100});
    auto cache = detail::init_cache<InputElementType, decltype(pred)>(cache_hint);
    std::for_each(first, last, [&] (const auto& value) {
        visitor(value, detail::any_of(adaptor(value), pred, cache));
    });
}

template <typename InputIterator, typename UnaryPredicate, typename BinaryFunction>
void for_each_any_of(InputIterator first, InputIterator last,
                     UnaryPredicate&& pred,
                     BinaryFunction&& visitor)
{
    for_each_any_of(first, last, std::forward<UnaryPredicate>(pred), std::forward<BinaryFunction>(visitor), [] (const auto& g) { return g; });
}

template <typename InputIterator, typename UnaryPredicate, typename BinaryFunction, typename UnaryFunction>
void for_each_all_of(InputIterator first, InputIterator last,
                     UnaryPredicate&& pred,
                     BinaryFunction&& visitor,
                     UnaryFunction&& adaptor)
{
    using V = typename std::iterator_traits<InputIterator>::value_type;
    using InputGenotypeType = std::result_of_t<UnaryFunction(V)>;
    using InputElementType = typename InputGenotypeType::ElementType;
    static_assert(detail::is_haplotype_like_v<InputElementType>, "");
    if (first == last) return;
    const auto cache_hint = std::min(static_cast<std::size_t>(std::distance(first, last)), std::size_t {100});
    auto cache = detail::init_cache<InputElementType, decltype(pred)>(cache_hint);
    std::for_each(first, last, [&] (const auto& value) {
        visitor(value, detail::all_of(adaptor(value), pred, cache));
    });
}

template <typename InputIterator, typename UnaryPredicate, typename BinaryFunction>
void for_each_all_of(InputIterator first, InputIterator last,
                     UnaryPredicate&& pred,
                     BinaryFunction&& visitor)
{
    for_each_all_of(first, last, std::forward<UnaryPredicate>(pred), std::forward<BinaryFunction>(visitor), [] (const auto& g) { return g; });
}

template <typename InputIterator, typename MappableType, typename BinaryFunction, typename UnaryFunction>
void for_each_contains(InputIterator first, InputIterator last,
                       const MappableType& target,
                       BinaryFunction&& visitor,
                       UnaryFunction&& adaptor)
{
    static_assert(!detail::is_haplotype_like_v<MappableType>, "");
    return for_each_any_of(first, last, [&] (const auto& element) { return contains(element, target); },
                           std::forward<BinaryFunction>(visitor), std::forward<UnaryFunction>(adaptor));
}

template <typename InputIterator, typename MappableType, typename BinaryFunction>
void for_each_contains(InputIterator first, InputIterator last,
                       const MappableType& target,
                       BinaryFunction&& visitor)
{
    for_each_contains(first, last, target, std::forward<BinaryFunction>(visitor), [] (const auto& g) { return g; });
}

template <typename InputIterator, typename MappableType, typename OutputIterator, typename BinaryFunction, typename UnaryFunction>
void transform_contains(InputIterator first, InputIterator last,
                        const MappableType& target,
                        OutputIterator result,
                        BinaryFunction&& f,
                        UnaryFunction&& adaptor)
{
    for_each_contains(first, last, target,
                      [&] (auto&& value, bool c) { *result++ = f(std::forward<decltype(value)>(value), c); },
                      std::forward<UnaryFunction>(adaptor));
}

namespace detail {

template <typename MappableType1, typename MappableType2>
bool contains(const Genotype<MappableType1>& lhs, const Genotype<MappableType2>& rhs, std::true_type)
{
    using std::cbegin; using std::cend; using std::begin; using std::end;
    using AlleleReference = std::reference_wrapper<const MappableType2>;
    if (lhs.ploidy() != rhs.ploidy()) return false;
    const auto lhs_copy = copy<MappableType2>(lhs, mapped_region(rhs));
    // Try to avoid sorting if possible
    if (std::is_sorted(cbegin(lhs_copy), cend(lhs_copy))) {
        if (std::is_sorted(cbegin(rhs), cend(rhs))) {
            return std::equal(cbegin(rhs), cend(rhs), cbegin(lhs_copy));
        }
        std::vector<AlleleReference> rhs_alleles {cbegin(rhs), cend(rhs)};
        std::sort(begin(rhs_alleles), end(rhs_alleles),
                  [] (const auto& lhs, const auto& rhs) { return lhs.get() < rhs.get(); });
        return std::equal(cbegin(lhs_copy), cend(lhs_copy), cbegin(rhs_alleles),
                          [] (const auto& lhs, const auto& rhs) { return lhs == rhs.get(); });
    }
    std::vector<AlleleReference> lhs_alleles {cbegin(lhs_copy), cend(lhs_copy)};
    std::sort(begin(lhs_alleles), end(lhs_alleles),
              [] (const auto& lhs, const auto& rhs) { return lhs.get() < rhs.get(); });
    
    if (std::is_sorted(cbegin(rhs), cend(rhs))) {
        return std::equal(cbegin(rhs), cend(rhs), cbegin(lhs_alleles),
                          [] (const auto& lhs, const auto& rhs) { return lhs == rhs.get(); });
    }
    std::vector<AlleleReference> rhs_alleles {cbegin(rhs), cend(rhs)};
    std::sort(begin(rhs_alleles), end(rhs_alleles),
              [] (const auto& lhs, const auto& rhs) { return lhs.get() < rhs.get(); });
    return std::equal(cbegin(lhs_alleles), cend(lhs_alleles), cbegin(rhs_alleles),
                      [] (const auto& lhs, const auto& rhs) { return lhs.get() == rhs.get(); });
}

template <typename MappableType1, typename MappableType2>
bool contains(const Genotype<MappableType1>& lhs, const Genotype<MappableType2>& rhs, std::false_type)
{
    return copy<MappableType2>(lhs, mapped_region(rhs)) == rhs;
}

} // namespace detail

template <typename MappableType1, typename MappableType2>
bool contains(const Genotype<MappableType1>& lhs, const Genotype<MappableType2>& rhs)
{
    using B = std::integral_constant<bool, detail::is_haplotype_like_v<MappableType1>
                                       && !detail::is_haplotype_like_v<MappableType2>>;
    return detail::contains(lhs, rhs, B {});
}

template <typename InputIterator, typename MappableType, typename BinaryFunction, typename UnaryFunction>
void for_each_contains(InputIterator first, InputIterator last,
                       const Genotype<MappableType>& target,
                       BinaryFunction&& visitor,
                       UnaryFunction&& adaptor)
{
    using V = typename std::iterator_traits<InputIterator>::value_type;
    using InputGenotypeType = std::result_of_t<UnaryFunction(V)>;
    using InputElementType = typename InputGenotypeType::ElementType;
    static_assert(detail::is_haplotype_like_v<InputElementType>, "");
    static_assert(!detail::is_haplotype_like_v<MappableType>, "");
    using MappableTypeRef = std::reference_wrapper<const MappableType>;
    if (first == last) return;
    std::vector<MappableTypeRef> sorted_target {std::cbegin(target), std::cend(target)};
    const static auto ref_less = [] (const auto& lhs, const auto& rhs) { return lhs.get() < rhs.get(); };
    std::sort(std::begin(sorted_target), std::end(sorted_target), ref_less);
    const auto copy_helper = [&] (const auto& element) { return copy<MappableType>(element, mapped_region(target)); };
    const auto cache_hint = std::min(static_cast<std::size_t>(std::distance(first, last) / ploidy(target)), std::size_t {100});
    auto cache = detail::init_cache<InputElementType, decltype(copy_helper)>(cache_hint);
    std::vector<MappableTypeRef> sorted_given {};
    std::for_each(first, last, [&] (const auto& value) {
        const auto& genotype = adaptor(value);
        bool does_contain {false};
        if (ploidy(genotype) == ploidy(target)) {
            const auto target_copy = detail::transform(genotype, copy_helper, cache);
            sorted_given.assign(std::cbegin(target_copy), std::cend(target_copy));
            std::sort(std::begin(sorted_given), std::end(sorted_given), ref_less);
            const static auto ref_equal = [] (const auto& lhs, const auto& rhs) { return lhs.get() == rhs.get(); };
            does_contain = std::equal(std::cbegin(sorted_given), std::cend(sorted_given), std::cbegin(sorted_target), ref_equal);
        }
        visitor(value, does_contain);
    });
}

template <typename InputIterator, typename MappableType, typename BinaryFunction>
void for_each_contains(InputIterator first, InputIterator last,
                       const Genotype<MappableType>& target,
                       BinaryFunction&& visitor)
{
    for_each_contains(first, last, target, std::forward<BinaryFunction>(visitor), [] (const auto& g) { return g; });
}

template <typename InputIterator, typename MappableType, typename OutputIterator, typename BinaryFunction, typename UnaryFunction>
void transform_contains(InputIterator first, InputIterator last,
                        const Genotype<MappableType>& target,
                        OutputIterator result,
                        BinaryFunction&& f,
                        UnaryFunction&& adaptor)
{
    for_each_contains(first, last, target,
                      [&] (auto&& value, bool c) { *result++ = f(std::forward<decltype(value)>(value), c); },
                      std::forward<UnaryFunction>(adaptor));
}

template <typename MappableType2, typename MappableType1>
bool are_equal_in_region(const Genotype<MappableType1>& lhs, const Genotype<MappableType1>& rhs,
                         const GenomicRegion& region)
{
    return copy<MappableType2>(lhs, region) == copy<MappableType2>(rhs, region);
}

template <typename InputIterator, typename UnaryPredicate, typename UnaryFunction>
bool all_of_all_of(InputIterator first, InputIterator last,
                   UnaryPredicate&& pred,
                   UnaryFunction&& adaptor)
{
    using V = typename std::iterator_traits<InputIterator>::value_type;
    using InputGenotypeType = std::result_of_t<UnaryFunction(V)>;
    using InputElementType = typename InputGenotypeType::ElementType;
    static_assert(detail::is_haplotype_like_v<InputElementType>, "");
    if (first == last) return false;
    const auto cache_hint = std::min(static_cast<std::size_t>(std::distance(first, last)), std::size_t {100});
    auto cache = detail::init_cache<InputElementType, decltype(pred)>(cache_hint);
    return std::all_of(first, last, [&] (const auto& value) { return detail::all_of(adaptor(value), pred, cache); });
}

template <typename InputIterator, typename UnaryPredicate>
bool all_of_all_of(InputIterator first, InputIterator last,
                   UnaryPredicate&& pred)
{
    return all_of_all_of(first, last, std::forward<UnaryPredicate>(pred), [] (const auto& g) { return g; });
}

template <typename MappableType>
bool is_homozygous(const Genotype<MappableType>& genotype, const MappableType& element)
{
    return count(genotype, element) == ploidy(genotype);
}
template <typename MappableType, std::enable_if_t<detail::is_haplotype_like_v<MappableType>, int> = 0>
bool is_homozygous(const Genotype<MappableType>& genotype, const GenomicRegion& region)
{
    return is_homozygous(copy<Allele>(genotype, region));
}
template <typename MappableType, std::enable_if_t<detail::is_haplotype_like_v<MappableType>, int> = 0>
bool is_homozygous(const Genotype<MappableType>& genotype, const Allele& allele)
{
    return all_of(genotype, [&] (const auto& haplotype) { return contains(haplotype, allele); });
}
template <typename MappableType, std::enable_if_t<detail::is_haplotype_like_v<MappableType>, int> = 0>
bool is_heterozygous(const Genotype<MappableType>& genotype, const GenomicRegion& region)
{
    return is_heterozygous(copy<Allele>(genotype, region));
}
template <typename MappableType, std::enable_if_t<detail::is_haplotype_like_v<MappableType>, int> = 0>
bool is_heterozygous(const Genotype<MappableType>& genotype, const Allele& allele)
{
    auto occ = count_if(genotype, [&] (const auto& haplotype) { return contains(haplotype, allele); });
    return occ > 0 && occ < ploidy(genotype);
}

template <typename InputIterator, typename MappableType, typename BinaryFunction, typename UnaryFunction>
void for_each_is_homozygous(InputIterator first, InputIterator last,
                            const MappableType& target,
                            BinaryFunction&& visitor,
                            UnaryFunction&& adaptor)
{
    for_each_all_of(first, last, [&] (const auto& element) { return contains(element, target); },
                    std::forward<BinaryFunction>(visitor), std::forward<UnaryFunction>(adaptor));
}
template <typename InputIterator, typename MappableType, typename BinaryFunction>
void for_each_is_homozygous(InputIterator first, InputIterator last,
                            const MappableType& target,
                            BinaryFunction&& visitor)
{
    for_each_homozygous(first, last, target, std::forward<BinaryFunction>(visitor), [] (const auto& g) { return g; });
}

template <typename InputIterator, typename MappableType, typename OutputIterator, typename BinaryFunction, typename UnaryFunction>
void transform_is_homozygous(InputIterator first, InputIterator last,
                             const MappableType& target,
                             OutputIterator result,
                             BinaryFunction&& f,
                             UnaryFunction&& adaptor)
{
    for_each_is_homozygous(first, last, target,
                          [&] (auto&& value, bool c) { *result++ = f(std::forward<decltype(value)>(value), c); },
                          std::forward<UnaryFunction>(adaptor));
}

template <typename InputIterator, typename UnaryFunction>
bool all_of_is_homozygous(InputIterator first, InputIterator last,
                          const Allele& target,
                          UnaryFunction&& adaptor)
{
    return all_of_all_of(first, last, [&] (const auto& element) { return contains(element, target); },
                         std::forward<UnaryFunction>(adaptor));
}

Genotype<Haplotype> remap(const Genotype<Haplotype>& genotype, const GenomicRegion& region);

template <typename MappableType>
bool operator==(const Genotype<MappableType>& lhs, const Genotype<MappableType>& rhs)
{
    return lhs.ploidy() == rhs.ploidy() && std::equal(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs));
}

struct GenotypeLess
{
    template <typename T>
    bool operator()(const Genotype<T>& lhs, const Genotype<T>& rhs) const
    {
        return std::lexicographical_compare(std::cbegin(lhs), std::cend(lhs),
                                            std::cbegin(rhs), std::cend(rhs));
    }
};

struct GenotypeHash
{
    template <typename T>
    std::size_t operator()(const Genotype<T>& genotype) const
    {
        return boost::hash_range(std::cbegin(genotype), std::cend(genotype));
    }
};

std::size_t num_genotypes(unsigned num_elements, unsigned ploidy);
boost::optional<std::size_t> num_genotypes_noexcept(unsigned num_elements, unsigned ploidy) noexcept;
std::size_t max_num_elements(std::size_t num_genotypes, unsigned ploidy);
std::size_t element_cardinality_in_genotypes(unsigned num_elements, unsigned ploidy);

template <typename MappableType>
unsigned count_shared(const Genotype<MappableType>& lhs, const Genotype<MappableType>& rhs)
{
    if (lhs.ploidy() <= rhs.ploidy()) {
        if (lhs.ploidy() < 2 || (lhs.ploidy() == 2 && is_heterozygous(lhs))) {
            return std::count_if(std::cbegin(lhs), std::cend(lhs),
                                 [&rhs] (const auto& element) { return contains(rhs, element); });
        } else {
            const auto lhs_unique = collapse(lhs);
            return std::count_if(std::cbegin(lhs_unique), std::cend(lhs_unique),
                                 [&rhs] (const auto& element) { return contains(rhs, element); });
        }
    } else {
        return count_shared(rhs, lhs);
    }
}

template <typename MappableType>
bool have_shared(const Genotype<MappableType>& lhs, const Genotype<MappableType>& rhs)
{
    if (lhs.ploidy() <= rhs.ploidy()) {
        return std::any_of(std::cbegin(lhs), std::cend(lhs), [&rhs] (const auto& element) { return contains(rhs, element); });
    } else {
        return have_shared(rhs, lhs);
    }
}

namespace detail {

template <typename T>
struct ValueType
{
    using type = std::decay_t<T>;
};
template <typename T>
struct ValueType<std::shared_ptr<T>>
{
    using type = std::decay_t<T>;
};
template <typename T>
struct ValueType<std::reference_wrapper<T>>
{
    using type = std::decay_t<T>;
};

template <typename T>
using value_type_t = typename ValueType<T>::type;

template <typename Range>
using GenotypeType = Genotype<value_type_t<typename Range::value_type>>;

template <typename Range>
auto construct_empty_genotype_container(const Range& elements)
{
    return std::vector<GenotypeType<Range>> {};
}
template <typename MappableType>
auto construct_empty_genotype_container(const MappableBlock<MappableType>& elements)
{
    return MappableBlock<Genotype<MappableType>> {mapped_region(elements)};
}

template <typename Range>
using ResultContainerType = std::remove_cv_t<decltype(detail::construct_empty_genotype_container(std::declval<Range>()))>;

template <typename Range>
auto generate_all_haploid_genotypes(const Range& elements)
{
    auto result = construct_empty_genotype_container(elements);
    result.reserve(elements.size());
    for (const auto& element : elements) {
        result.emplace_back(1, element);
    }
    return result;
}

template <typename Range>
auto generate_all_diploid_biallelic_genotypes(const Range& elements)
{
    using GenotypeTp = GenotypeType<Range>;
    using ResultType = decltype(construct_empty_genotype_container(elements));
    return ResultType {
    GenotypeTp {elements[0], elements[0]},
    GenotypeTp {elements[0], elements[1]},
    GenotypeTp {elements[1], elements[1]}
    };
}

template <typename Range>
auto generate_all_diploid_triallelic_genotypes(const Range& elements)
{
    using GenotypeTp = GenotypeType<Range>;
    using ResultType = decltype(construct_empty_genotype_container(elements));
    return ResultType {
    GenotypeTp {elements[0], elements[0]},
    GenotypeTp {elements[0], elements[1]},
    GenotypeTp {elements[0], elements[2]},
    GenotypeTp {elements[1], elements[1]},
    GenotypeTp {elements[1], elements[2]},
    GenotypeTp {elements[2], elements[2]}
    };
}

template <typename Range>
auto generate_all_diploid_tetraallelic_genotypes(const Range& elements)
{
    using GenotypeTp = GenotypeType<Range>;
    using ResultType = decltype(construct_empty_genotype_container(elements));
    return ResultType {
    GenotypeTp {elements[0], elements[0]},
    GenotypeTp {elements[0], elements[1]},
    GenotypeTp {elements[0], elements[2]},
    GenotypeTp {elements[0], elements[3]},
    GenotypeTp {elements[1], elements[1]},
    GenotypeTp {elements[1], elements[2]},
    GenotypeTp {elements[1], elements[3]},
    GenotypeTp {elements[2], elements[2]},
    GenotypeTp {elements[2], elements[3]},
    GenotypeTp {elements[3], elements[3]}
    };
}

template <typename Range>
auto generate_all_triploid_biallelic_genotypes(const Range& elements)
{
    using GenotypeTp = GenotypeType<Range>;
    using ResultType = decltype(construct_empty_genotype_container(elements));
    return ResultType {
    GenotypeTp {elements[0], elements[0], elements[0]},
    GenotypeTp {elements[0], elements[0], elements[1]},
    GenotypeTp {elements[0], elements[1], elements[1]},
    GenotypeTp {elements[1], elements[1], elements[1]}
    };
}

template <typename Range>
auto generate_genotype(const Range& elements, const std::vector<unsigned>& element_indicies)
{
    GenotypeType<Range> result {static_cast<unsigned>(element_indicies.size())};
    for (const auto i : element_indicies) {
        result.emplace(elements[i]);
    }
    return result;
}

} // namespace detail

template <typename Range>
auto generate_all_genotypes(const Range& elements, const unsigned ploidy)
{
    using GenotypeTp = detail::GenotypeType<Range>;
    using ResultType = detail::ResultContainerType<Range>;
    // Optimise simple cases
    if (ploidy == 0 || elements.empty()) {
        return detail::construct_empty_genotype_container(elements);
    }
    const auto num_elements = static_cast<unsigned>(elements.size());
    if (num_elements == 1) {
        return ResultType {GenotypeTp {ploidy, elements.front()}};
    }
    if (ploidy == 2) {
        switch (num_elements) {
            case 2:
                return detail::generate_all_diploid_biallelic_genotypes(elements);
            case 3:
                return detail::generate_all_diploid_triallelic_genotypes(elements);
            case 4:
                return detail::generate_all_diploid_tetraallelic_genotypes(elements);
        }
    }
    if (ploidy == 1) {
        return detail::generate_all_haploid_genotypes(elements);
    }
    if (ploidy == 3 && num_elements == 2) {
        return detail::generate_all_triploid_biallelic_genotypes(elements);
    }
    // Otherwise resort to general algorithm
    ResultType result {detail::construct_empty_genotype_container(elements)};
    result.reserve(num_genotypes(num_elements, ploidy));
    std::vector<unsigned> element_indicies(ploidy, 0);
    while (true) {
        if (element_indicies[0] == num_elements) {
            unsigned i {0};
            while (++i < ploidy && element_indicies[i] == num_elements - 1);
            if (i == ploidy) break;
            ++element_indicies[i];
            std::fill_n(std::begin(element_indicies), i + 1, element_indicies[i]);
        }
        result.push_back(detail::generate_genotype(elements, element_indicies));
        ++element_indicies[0];
    }
    return result;
}

template <typename Range, typename UnaryPredicate, typename OutputIterator>
auto generate_all_genotypes(const Range& elements, const unsigned ploidy,
                            UnaryPredicate&& pred, OutputIterator result_itr)
{
    if (ploidy == 0 || elements.empty()) return result_itr;
    const auto num_elements = static_cast<unsigned>(elements.size());
    std::vector<unsigned> element_indicies(ploidy, 0);
    while (true) {
        if (element_indicies[0] == num_elements) {
            unsigned i {0};
            while (++i < ploidy && element_indicies[i] == num_elements - 1);
            if (i == ploidy) break;
            ++element_indicies[i];
            std::fill_n(std::begin(element_indicies), i + 1, element_indicies[i]);
        }
        auto genotype = detail::generate_genotype(elements, element_indicies);
        if (pred(genotype)) *result_itr++ = std::move(genotype);
        ++element_indicies[0];
    }
    return result_itr;
}

std::size_t num_max_zygosity_genotypes(unsigned num_elements, unsigned ploidy);
boost::optional<std::size_t> num_max_zygosity_genotypes_noexcept(unsigned num_elements, unsigned ploidy) noexcept;

template <typename Range, typename OutputIterator>
void
generate_all_max_zygosity_genotypes(const Range& elements, const unsigned ploidy, OutputIterator result)
{
    generate_all_genotypes(elements, ploidy, [] (const auto& genotype) { return is_max_zygosity(genotype); }, result);
}

template <typename Range>
auto generate_all_max_zygosity_genotypes(const Range& elements, const unsigned ploidy)
{
    auto result = detail::construct_empty_genotype_container(elements);
    if (elements.size() >= ploidy) {
        try {
            result.reserve(num_max_zygosity_genotypes(elements.size(), ploidy));
        } catch (const std::overflow_error& e) {
            result.reserve(std::numeric_limits<std::size_t>::max()); // Probably this will throw a bad_alloc
        }
    }
    generate_all_max_zygosity_genotypes(elements, ploidy, std::back_inserter(result));
    return result;
}

template <typename ForwardIterator, typename MappableType, 
          typename OutputIterator, typename UnaryPredicate>
OutputIterator
extend(ForwardIterator first, ForwardIterator last, 
       const MappableBlock<MappableType>& haplotypes,
       OutputIterator result,
       UnaryPredicate&& selector)
{
    std::for_each(first, last, [&] (const auto& genotype) {
        for (const auto& haplotype : haplotypes) {
            if (selector(genotype, haplotype)) {
                auto extended_genotype = genotype;
                extended_genotype.emplace(haplotype);
                extended_genotype.shrink_to_fit();
                *result++ = std::move(extended_genotype);
            }
        }
    });
    return result;
}

template <typename ForwardIterator, typename MappableType, 
          typename OutputIterator>
OutputIterator
extend(ForwardIterator first, ForwardIterator last, 
       const MappableBlock<MappableType>& haplotypes,
       OutputIterator result)
{
    const static auto default_selector = [] (const auto&, const auto&) noexcept { return true; };
    return extend(first, last, haplotypes, result, default_selector);
}

template <typename MappableType, typename UnaryPredicate>
auto extend(const MappableBlock<Genotype<MappableType>>& genotypes,
            const MappableBlock<MappableType>& haplotypes,
            UnaryPredicate&& selector)
{
    MappableBlock<Genotype<MappableType>> result {mapped_region(haplotypes)};
    result.reserve(genotypes.size() * haplotypes.size());
    extend(std::cbegin(genotypes), std::cend(genotypes), haplotypes, std::back_inserter(result), std::forward<UnaryPredicate>(selector));
    return result;
}

template <typename MappableType>
auto extend(const MappableBlock<Genotype<MappableType>>& genotypes,
            const MappableBlock<MappableType>& haplotypes)
{
    const static auto default_selector = [] (const auto&, const auto&) noexcept { return true; };
    return extend(genotypes, haplotypes, default_selector);
}

namespace detail {

inline std::size_t estimate_num_elements(const std::size_t num_genotypes)
{
    return num_genotypes;
}

} // namespace detail

template <typename Range>
auto extract_unique_elements(const Range& genotypes)
{
    using MappableType = typename Range::value_type::ElementType;
    std::unordered_set<std::reference_wrapper<const MappableType>> unique_elements {};
    unique_elements.reserve(detail::estimate_num_elements(genotypes.size()));
    for (const auto& genotype : genotypes) {
        for (const auto& element : genotype) {
            unique_elements.emplace(element);
        }
    }
    std::vector<MappableType> result {};
    result.reserve(unique_elements.size());
    std::transform(std::cbegin(unique_elements), std::cend(unique_elements), std::back_inserter(result),
                   [] (const auto& element_ref) { return element_ref.get(); });
    return result;
}

template <typename Range>
auto extract_unique_element_refs(const Range& genotypes)
{
    using MappableType = typename Range::value_type::ElementType;
    std::unordered_set<std::reference_wrapper<const MappableType>> unique_elements {};
    unique_elements.reserve(genotypes.size());
    for (const auto& genotype : genotypes) {
        for (const auto& element : genotype) {
            unique_elements.emplace(element);
        }
    }
    std::vector<std::reference_wrapper<const MappableType>> result {};
    result.reserve(unique_elements.size());
    std::copy(std::cbegin(unique_elements), std::cend(unique_elements), std::back_inserter(result));
    return result;
}

template <typename MappableType>
auto make_element_count_map(const Genotype<MappableType>& genotype)
{
    std::unordered_map<MappableType, unsigned> result {};
    result.reserve(zygosity(genotype));
    for (unsigned i {0}; i < genotype.ploidy(); ++i) {
        ++result[genotype[i]];
    }
    return result;
}

template <typename MappableType>
auto make_element_ref_count_map(const Genotype<MappableType>& genotype)
{
    std::unordered_map<std::reference_wrapper<const MappableType>, unsigned> result {};
    result.reserve(genotype.zygosity());
    for (unsigned i {0}; i < genotype.ploidy(); ++i) {
        ++result[genotype.at(i)];
    }
    return result;
}

template <typename MappableType2, typename Container,
          typename = std::enable_if_t<!std::is_same<typename Container::value_type, Haplotype>::value>>
auto copy_unique(const Container& genotypes, const GenomicRegion& region)
{
    std::vector<Genotype<MappableType2>> result {};
    result.reserve(genotypes.size());
    for (const auto& genotype : genotypes) {
        result.push_back(copy<MappableType2>(genotype, region));
    }
    std::sort(std::begin(result), std::end(result), GenotypeLess {});
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

template <typename MappableType>
std::ostream& operator<<(std::ostream& os, const Genotype<MappableType>& genotype)
{
    if (genotype.ploidy() == 0) {
        os << "empty genotype";
        return os;
    }
    const auto element_counts = make_element_count_map(genotype);
    std::vector<std::pair<MappableType, unsigned>> p {element_counts.begin(), element_counts.end()};
    for (unsigned i {0}; i < p.size() - 1; ++i) {
        os << p[i].first << "(" << p[i].second << "),";
    }
    os << p.back().first << "(" << p.back().second << ")";
    return os;
}

namespace debug {

template <typename S, typename MappableType>
void print_alleles(S&& stream, const Genotype<MappableType>& genotype)
{
    if (genotype.ploidy() == 0) {
        stream << "[]";
    }
    const auto haplotype_counts = make_element_count_map(genotype);
    std::vector<std::pair<MappableType, unsigned>> p {haplotype_counts.begin(), haplotype_counts.end()};
    stream << "[";
    for (unsigned i {0}; i < p.size(); ++i) {
        print_alleles(stream, p[i].first);
        stream << "(" << p[i].second << ")";
        if (i < p.size() - 1) stream << ",";
    }
    stream << "]";
}

template <typename S, typename MappableType>
void print_variant_alleles(S&& stream, const Genotype<MappableType>& genotype)
{
    if (genotype.ploidy() == 0) {
        stream << "[]";
    }
    const auto simplified = collapse(genotype);
    stream << "[";
    for (unsigned i {0}; i < simplified.ploidy(); ++i) {
        print_variant_alleles(stream, simplified[i]);
        stream << "(" << count(genotype, simplified[i]) << ")";
        if (i < simplified.ploidy() - 1) stream << ",";
    }
    stream << "]";
}

template <typename MappableType>
void print_alleles(const Genotype<MappableType>& genotype)
{
    print_alleles(std::cout, genotype);
}
template <typename MappableType>
void print_variant_alleles(const Genotype<MappableType>& genotype)
{
    print_variant_alleles(std::cout, genotype);
}

} // namespace debug

} // namespace octopus

namespace std {

template <typename MappableType> struct hash<octopus::Genotype<MappableType>>
{
    size_t operator()(const octopus::Genotype<MappableType>& genotype) const
    {
        return octopus::GenotypeHash()(genotype);
    }
};

template <typename MappableType> struct hash<reference_wrapper<const octopus::Genotype<MappableType>>>
{
    size_t operator()(const reference_wrapper<const octopus::Genotype<MappableType>> genotype) const
    {
        return hash<octopus::Genotype<MappableType>>()(genotype);
    }
};

} // namespace std

#endif
