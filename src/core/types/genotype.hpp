// Copyright (c) 2015-2019 Daniel Cooke
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

#include <boost/functional/hash.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "concepts/equitable.hpp"
#include "concepts/mappable.hpp"
#include "containers/mappable_block.hpp"
#include "utils/reorder.hpp"
#include "allele.hpp"
#include "haplotype.hpp"
#include "indexed_haplotype.hpp"

namespace octopus {

namespace detail {

template <typename T>
constexpr bool is_allele = std::is_same<T, Allele>::value || std::is_same<T, ContigAllele>::value;
template <typename T>
constexpr bool is_haplotype_or_indexed_haplotype = std::is_same<T, Haplotype>::value || is_indexed_haplotype<T>;

} // namespace detail

template <typename T>
constexpr bool is_genotypeable = detail::is_allele<T> || detail::is_haplotype_or_indexed_haplotype<T>;

//template <typename MappableType, typename = std::enable_if_t<is_genotypeable<MappableType>>> class Genotype;

template <typename MappableType, typename Enable = void>
class Genotype : public Equitable<Genotype<MappableType>>, public Mappable<Genotype<MappableType>>
{
public:
    static_assert(is_genotypeable<MappableType>, "");
    
    using ElementType   = MappableType;
    using MappingDomain = RegionType<ElementType>;
    
    using ordered = std::false_type;
    using share_memory = std::false_type;
    
    Genotype() = default;
    
    explicit Genotype(unsigned ploidy);
    explicit Genotype(unsigned ploidy, const MappableType& init);
    explicit Genotype(std::initializer_list<MappableType> alleles);
    
    Genotype(const Genotype&)            = default;
    Genotype& operator=(const Genotype&) = default;
    Genotype(Genotype&&)                 = default;
    Genotype& operator=(Genotype&&)      = default;
    
    ~Genotype() = default;
    
    template <typename T> void emplace(T&& element);
    
    const MappingDomain& mapped_region() const noexcept;
    
    const MappableType& operator[](unsigned n) const;
    
    unsigned ploidy() const noexcept;
    
    bool contains(const MappableType& element) const;
    unsigned count(const MappableType& element) const;
    
    bool is_homozygous() const;
    unsigned zygosity() const;
    
    Genotype<MappableType> collapse() const;
    
    std::vector<MappableType> copy_unique() const;
    
    void reorder_alleles(const std::vector<unsigned>& order);
    
private:
    std::vector<MappableType> elements_;

public:
    using Iterator = typename decltype(elements_)::const_iterator;
    
    Iterator begin() const noexcept { return std::begin(elements_); }
    Iterator end() const noexcept { return std::end(elements_); }
    Iterator cbegin() const noexcept { return std::cbegin(elements_); }
    Iterator cend() const noexcept { return std::cend(elements_); }
};

namespace detail {

template <typename MappableType>
using enable_if_haplotypelike = std::enable_if_t<is_haplotype_or_indexed_haplotype<MappableType>>;

} // namespace detail

template <typename MappableType>
class Genotype<MappableType, detail::enable_if_haplotypelike<MappableType>>
    : public Equitable<Genotype<MappableType>>, public Mappable<Genotype<MappableType>>
{
public:
    using ElementType   = MappableType;
    using MappingDomain = typename MappableType::MappingDomain;
    
    using ordered = std::true_type;
    using share_memory = std::true_type;
    
    Genotype() = default;
    
    explicit Genotype(unsigned ploidy);
    explicit Genotype(unsigned ploidy, const ElementType& init);
    explicit Genotype(unsigned ploidy, const std::shared_ptr<ElementType>& init);
    explicit Genotype(std::initializer_list<ElementType> elements);
    explicit Genotype(std::initializer_list<std::shared_ptr<ElementType>> elements);
    
    Genotype(const Genotype&)            = default;
    Genotype& operator=(const Genotype&) = default;
    Genotype(Genotype&&)                 = default;
    Genotype& operator=(Genotype&&)      = default;
    
    ~Genotype() = default;
    
    template <typename T> void emplace(T&& element);
    void emplace(const std::shared_ptr<ElementType>& element);
    
    const ElementType& operator[](unsigned n) const;
    
    const GenomicRegion& mapped_region() const noexcept;
    
    unsigned ploidy() const noexcept;
    
    bool contains(const ElementType& haplotype) const;
    unsigned count(const ElementType& haplotype) const;
    
    bool is_homozygous() const;
    unsigned zygosity() const;
    
    Genotype<ElementType> collapse() const;
    
    std::vector<ElementType> copy_unique() const;
    std::vector<std::reference_wrapper<const ElementType>> copy_unique_ref() const;
    std::vector<unsigned> unique_counts() const;
    
private:
    using HaplotypePtr  = std::shared_ptr<ElementType>;
    using BaseContainer = std::vector<HaplotypePtr>;
    using BaseIterator  = typename BaseContainer::const_iterator;
    
    BaseContainer haplotypes_;
    
    struct HaplotypePtrLess
    {
        bool operator()(const HaplotypePtr& lhs, const HaplotypePtr& rhs) const { return *lhs < *rhs; }
        bool operator()(const ElementType& lhs, const HaplotypePtr& rhs) const { return lhs < *rhs; }
        bool operator()(const HaplotypePtr& lhs, const ElementType& rhs) const  { return *lhs < rhs; }
    };
    
    struct HaplotypePtrEqual
    {
        bool operator()(const HaplotypePtr& lhs, const HaplotypePtr& rhs) const { return lhs == rhs || *lhs == *rhs; }
        bool operator()(const ElementType& lhs, const HaplotypePtr& rhs) const { return lhs == *rhs; }
        bool operator()(const HaplotypePtr& lhs, const ElementType& rhs) const { return *lhs == rhs; }
    };
    
public:
    class Iterator : public BaseIterator
    {
    public:
        using value_type = typename BaseIterator::value_type::element_type;
        using reference  = const value_type&;
        using pointer    = const value_type*;
        
        Iterator(BaseIterator it) : BaseIterator {it} {}
        reference operator*() const noexcept { return *BaseIterator::operator*(); }
        pointer operator->() const noexcept { return std::addressof(this->operator*()); }
    };
    
    Iterator begin() const noexcept { return std::begin(haplotypes_); }
    Iterator end() const noexcept { return std::end(haplotypes_); }
    Iterator cbegin() const noexcept { return std::cbegin(haplotypes_); }
    Iterator cend() const noexcept { return std::cend(haplotypes_); }
};

// Genotype<MappableType>

template <typename MappableType, typename Enable>
template <typename T>
void Genotype<MappableType, Enable>::emplace(T&& element)
{
    elements_.emplace_back(std::forward<T>(element));
}

template <typename MappableType, typename Enable>
Genotype<MappableType, Enable>::Genotype(const unsigned ploidy)
: elements_ {}
{
    elements_.reserve(ploidy);
}

template <typename MappableType, typename Enable>
Genotype<MappableType, Enable>::Genotype(const unsigned ploidy, const MappableType& init)
: elements_ {ploidy, init}
{}

template <typename MappableType, typename Enable>
Genotype<MappableType, Enable>::Genotype(std::initializer_list<MappableType> elements)
: elements_ {elements}
{}

template <typename MappableType, typename Enable>
const typename Genotype<MappableType, Enable>::MappingDomain& Genotype<MappableType, Enable>::mapped_region() const noexcept
{
    return elements_.front().mapped_region();
}

template <typename MappableType, typename Enable>
const MappableType& Genotype<MappableType, Enable>::operator[](const unsigned n) const
{
    return elements_[n];
}

template <typename MappableType, typename Enable>
unsigned Genotype<MappableType, Enable>::ploidy() const noexcept
{
    return elements_.size();
}

template <typename MappableType, typename Enable>
bool Genotype<MappableType, Enable>::is_homozygous() const
{
    return std::adjacent_find(std::cbegin(elements_), std::cend(elements_),
                              std::not_equal_to<Allele>()) == std::cend(elements_);
}

template <typename MappableType, typename Enable>
unsigned Genotype<MappableType, Enable>::zygosity() const
{
    if (ploidy() == 1 || is_homozygous()) {
        return 1;
    } else if (ploidy() == 2) {
        return 2;
    }
    return copy_unique().size();
}

template <typename MappableType, typename Enable>
bool Genotype<MappableType, Enable>::contains(const MappableType& element) const
{
    return std::find(std::cbegin(elements_), std::cend(elements_), element) != std::cend(elements_);
}

template <typename MappableType, typename Enable>
unsigned Genotype<MappableType, Enable>::count(const MappableType& element) const
{
    return std::count(std::cbegin(elements_), std::cend(elements_), element);
}

template <typename MappableType, typename Enable>
Genotype<MappableType> Genotype<MappableType, Enable>::collapse() const
{
    if (ploidy() < 2) return *this;
    Genotype<MappableType> result {ploidy()};
    result.emplace(elements_.front());
    std::for_each(std::next(std::cbegin(elements_)), std::cend(elements_), [&] (const auto& element) {
        if (!result.contains(element)) {
            result.emplace(element);
        }
    });
    return result;
}

template <typename MappableType, typename Enable>
std::vector<MappableType> Genotype<MappableType, Enable>::copy_unique() const
{
    auto result = elements_;
    std::sort(std::begin(result), std::end(result));
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

template <typename MappableType, typename Enable>
void Genotype<MappableType, Enable>::reorder_alleles(const std::vector<unsigned>& order)
{
    assert(order.size() == elements_.size());
    reorder(std::cbegin(order), std::cend(order), std::begin(elements_));
}

// Genotype<IndexedHaplotype>

template <typename MappableType>
Genotype<MappableType, detail::enable_if_haplotypelike<MappableType>>::Genotype(const unsigned ploidy)
: haplotypes_ {}
{
    haplotypes_.reserve(ploidy);
}

template <typename MappableType>
Genotype<MappableType, detail::enable_if_haplotypelike<MappableType>>::Genotype(const unsigned ploidy, const MappableType& init)
{
    if (ploidy > 0) {
        haplotypes_.resize(ploidy);
        haplotypes_.front() = std::make_shared<Haplotype>(init);
        std::fill_n(std::next(std::begin(haplotypes_)), ploidy - 1, haplotypes_.front());
    }
}

template <typename MappableType>
Genotype<MappableType, detail::enable_if_haplotypelike<MappableType>>::Genotype(const unsigned ploidy, const std::shared_ptr<MappableType>& init)
: haplotypes_ {ploidy, init}
{}

template <typename MappableType>
Genotype<MappableType, detail::enable_if_haplotypelike<MappableType>>::Genotype(std::initializer_list<MappableType> haplotypes)
{
    std::transform(std::cbegin(haplotypes), std::cend(haplotypes), std::back_inserter(haplotypes_),
                   [] (const auto& haplotype) { return std::make_shared<Haplotype>(haplotype); });
    std::sort(std::begin(haplotypes_), std::end(haplotypes_), HaplotypePtrLess {});
}

template <typename MappableType>
Genotype<MappableType, detail::enable_if_haplotypelike<MappableType>>::Genotype(std::initializer_list<std::shared_ptr<MappableType>> haplotypes)
: haplotypes_ {haplotypes}
{
    std::sort(std::begin(haplotypes_), std::end(haplotypes_), HaplotypePtrLess {});
}

template <typename MappableType>
template <typename T>
void Genotype<MappableType, detail::enable_if_haplotypelike<MappableType>>::emplace(T&& haplotype)
{
    haplotypes_.emplace_back(std::make_shared<MappableType>(std::forward<T>(haplotype)));
    std::inplace_merge(std::begin(haplotypes_), std::prev(std::end(haplotypes_)),
                       std::end(haplotypes_), HaplotypePtrLess {});
}

namespace detail {

template <typename T, typename BinaryPredicate>
typename std::vector<T>::iterator
insert_sorted(T value, std::vector<T>& values, BinaryPredicate compare)
{
    auto position = std::upper_bound(std::begin(values), std::end(values), value, compare);
    return values.insert(position, std::move(value));
}

} // namespace detail

template <typename MappableType>
void Genotype<MappableType, detail::enable_if_haplotypelike<MappableType>>::emplace(const std::shared_ptr<MappableType>& haplotype)
{
    detail::insert_sorted(haplotype, haplotypes_, HaplotypePtrLess {});
}

template <typename MappableType>
const MappableType& Genotype<MappableType, detail::enable_if_haplotypelike<MappableType>>::operator[](const unsigned n) const
{
    return *haplotypes_[n];
}

template <typename MappableType>
const GenomicRegion& Genotype<MappableType, detail::enable_if_haplotypelike<MappableType>>::mapped_region() const noexcept
{
    return haplotypes_.front()->mapped_region();
}

template <typename MappableType>
unsigned Genotype<MappableType, detail::enable_if_haplotypelike<MappableType>>::ploidy() const noexcept
{
    return static_cast<unsigned>(haplotypes_.size());
}

template <typename MappableType>
bool Genotype<MappableType, detail::enable_if_haplotypelike<MappableType>>::is_homozygous() const
{
    return haplotypes_.front() == haplotypes_.back() ||  *haplotypes_.front() == *haplotypes_.back();
}

template <typename MappableType>
unsigned Genotype<MappableType, detail::enable_if_haplotypelike<MappableType>>::zygosity() const
{
    if (ploidy() < 2) return ploidy();
    unsigned result {0};
    for (auto it = std::cbegin(haplotypes_), last = std::cend(haplotypes_); it != last; ++result) {
        // naive algorithm faster in practice than binary searching
        it = std::find_if_not(std::next(it), last, [it] (const auto& x) { return x == *it || *x == **it; });
    }
    return result;
}

template <typename MappableType>
Genotype<MappableType>
Genotype<MappableType, detail::enable_if_haplotypelike<MappableType>>::collapse() const
{
    if (zygosity() == ploidy()) return *this;
    std::vector<std::reference_wrapper<const HaplotypePtr>> unique {};
    unique.reserve(ploidy());
    std::unique_copy(std::cbegin(haplotypes_), std::cend(haplotypes_), std::back_inserter(unique), HaplotypePtrEqual {});
    Genotype<Haplotype> result {static_cast<unsigned>(unique.size())};
    for (const auto& haplotype_ptr : unique) result.emplace(haplotype_ptr.get());
    return result;
}

template <typename MappableType>
bool Genotype<MappableType, detail::enable_if_haplotypelike<MappableType>>::contains(const MappableType& haplotype) const
{
    return std::binary_search(std::cbegin(haplotypes_), std::cend(haplotypes_), haplotype, HaplotypePtrLess {});
}

template <typename MappableType>
unsigned Genotype<MappableType, detail::enable_if_haplotypelike<MappableType>>::count(const MappableType& haplotype) const
{
    const auto equal_range = std::equal_range(std::cbegin(haplotypes_), std::cend(haplotypes_),
                                              haplotype, HaplotypePtrLess {});
    return static_cast<unsigned>(std::distance(equal_range.first, equal_range.second));
}

template <typename MappableType>
std::vector<MappableType>
Genotype<MappableType, detail::enable_if_haplotypelike<MappableType>>::copy_unique() const
{
    std::vector<std::reference_wrapper<const HaplotypePtr>> ptr_copy {};
    ptr_copy.reserve(ploidy());
    std::unique_copy(std::cbegin(haplotypes_), std::cend(haplotypes_), std::back_inserter(ptr_copy), HaplotypePtrEqual {});
    std::vector<Haplotype> result {};
    result.reserve(ptr_copy.size());
    std::transform(std::cbegin(ptr_copy), std::cend(ptr_copy), std::back_inserter(result),
                   [] (const auto& ptr) { return *ptr.get(); });
    return result;
}

template <typename MappableType>
std::vector<std::reference_wrapper<const MappableType>>
Genotype<MappableType, detail::enable_if_haplotypelike<MappableType>>::copy_unique_ref() const
{
    std::vector<std::reference_wrapper<const MappableType>> result {};
    result.reserve(ploidy());
    std::transform(std::cbegin(haplotypes_), std::cend(haplotypes_), std::back_inserter(result),
                   [] (const HaplotypePtr& haplotype) { return std::cref(*haplotype); });
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

template <typename MappableType>
std::vector<unsigned>
Genotype<MappableType, detail::enable_if_haplotypelike<MappableType>>::unique_counts() const
{
    std::vector<unsigned> result {};
    result.reserve(haplotypes_.size());
    for (auto itr = std::cbegin(haplotypes_), last = std::cend(haplotypes_); itr != last;) {
        auto next = std::find_if_not(std::next(itr), last, [itr] (const auto& x) { return *x == **itr; });
        result.push_back(std::distance(itr, next));
        itr = next;
    }
    return result;
}

// non-member methods

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

template <typename T>
decltype(auto) get(const Genotype<T>& genotype) noexcept
{
    return genotype;
}

template <typename T>
decltype(auto) get(const T& genotype) noexcept
{
    return genotype.get();
}

} // namespace detail

template <typename MappableType2, typename MappableType1>
Genotype<MappableType2> convert(const Genotype<MappableType1>& genotype)
{
    Genotype<MappableType2> result {genotype.ploidy()};
    for (const auto& mappable : genotype) {
        result.emplace(mappable);
    }
    return result;
}

template <typename MappableType, typename G>
Genotype<MappableType> copy(const G& genotype, const GenomicRegion& region)
{
    Genotype<MappableType> result {detail::get(genotype).ploidy()};
    for (const auto& element : detail::get(genotype)) {
        result.emplace(copy<MappableType>(element, region));
    }
    return result;
}

namespace detail {

template <typename MappableType, typename Container>
auto copy_each_basic(const Container& genotypes, const GenomicRegion& region)
{
    std::vector<Genotype<MappableType>> result {};
    result.reserve(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::back_inserter(result),
                   [&] (const auto& genotype) { return copy<MappableType>(genotype, region); });
    return result;
}

template <typename MappableType>
struct share_memory : public Genotype<MappableType>::share_memory {};

template <typename T>
struct CopyType
{
    using type = typename T::ElementType;
};

template <typename T>
struct CopyType<std::reference_wrapper<const T>>
{
    using type = typename T::ElementType;
};

template <typename T>
using copy_type = typename CopyType<T>::type;

template <typename MappableType, typename G, typename Map>
Genotype<MappableType> copy(const G& genotype, const GenomicRegion& region, Map& cache)
{
    Genotype<MappableType> result {get(genotype).ploidy()};
    for (const auto& element : get(genotype)) {
        const auto itr = cache.find(element);
        if (itr == std::cend(cache)) {
            auto chunk = copy<MappableType>(element, region);
            cache.emplace(std::piecewise_construct, std::forward_as_tuple(element), std::forward_as_tuple(chunk));
            result.emplace(std::move(chunk));
        } else {
            result.emplace(itr->second);
        }
    }
    return result;
}

template <typename MappableType, typename Container>
auto copy_each_cached(const Container& genotypes, const GenomicRegion& region, std::false_type)
{
    using MappableType2 = copy_type<typename Container::value_type>;
    std::vector<Genotype<MappableType>> result {};
    if (genotypes.empty()) return result;
    std::unordered_map<MappableType2, MappableType> cache {};
    cache.reserve(genotypes.size());
    result.reserve(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::back_inserter(result),
                   [&] (const auto& genotype) { return copy<MappableType>(genotype, region, cache); });
    return result;
}

template <typename T>
void emplace(const std::shared_ptr<T>& element, Genotype<T>& genotype)
{
    genotype.emplace(element);
}

template <typename MappableType, typename G, typename Map>
Genotype<MappableType> copy_shared(const G& genotype, const GenomicRegion& region, Map& cache)
{
    Genotype<MappableType> result {get(genotype).ploidy()};
    for (const auto& element : get(genotype)) {
        const auto itr = cache.find(element);
        if (itr == std::cend(cache)) {
            auto copy_ptr = std::make_shared<MappableType>(copy<MappableType>(element, region));
            cache.emplace(std::piecewise_construct, std::forward_as_tuple(element), std::forward_as_tuple(copy_ptr));
            emplace(copy_ptr, result);
        } else {
            emplace(itr->second, result);
        }
    }
    return result;
}

template <typename MappableType, typename Container>
auto copy_each_cached(const Container& genotypes, const GenomicRegion& region, std::true_type)
{
    using MappableType2 = copy_type<typename Container::value_type>;
    std::vector<Genotype<MappableType>> result {};
    if (genotypes.empty()) return result;
    std::unordered_map<MappableType2, std::shared_ptr<MappableType>> cache {};
    cache.reserve(genotypes.size());
    result.reserve(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::back_inserter(result),
                   [&] (const auto& genotype) { return copy_shared<MappableType>(genotype, region, cache); });
    return result;
}

template <typename MappableType, typename Container>
auto copy_each_cached(const Container& genotypes, const GenomicRegion& region)
{
    using MappableType2 = copy_type<typename Container::value_type>;
    return copy_each_cached<MappableType>(genotypes, region, share_memory<MappableType2> {});
}

} // namespace detail

template <typename GenotypeType>
constexpr bool is_ordered = GenotypeType::ordered::value;

template <typename MappableType, typename Container>
std::vector<Genotype<MappableType>> copy_each(const Container& genotypes, const GenomicRegion& region)
{
    if (genotypes.size() < 10) {
        return detail::copy_each_basic<MappableType>(genotypes, region);
    } else {
        return detail::copy_each_cached<MappableType>(genotypes, region);
    }
}

template <typename UnaryPredicate>
bool any_of(const Genotype<Haplotype>& genotype, UnaryPredicate&& pred)
{
    for (auto haplotype_itr = std::cbegin(genotype); haplotype_itr != std::cend(genotype);) {
        if (pred(*haplotype_itr)) return true;
        const auto is_same_haplotype = [haplotype_itr] (const Haplotype& x) noexcept {
            return std::addressof(x) == std::addressof(*haplotype_itr); };
        haplotype_itr = std::find_if_not(std::next(haplotype_itr), std::cend(genotype), is_same_haplotype);
    }
    return false;
}

template <typename UnaryPredicate>
bool all_of(const Genotype<Haplotype>& genotype, UnaryPredicate&& pred)
{
    for (auto haplotype_itr = std::cbegin(genotype); haplotype_itr != std::cend(genotype);) {
        if (!pred(*haplotype_itr)) return false;
        const auto is_same_haplotype = [haplotype_itr] (const Haplotype& x) noexcept {
            return std::addressof(x) == std::addressof(*haplotype_itr); };
        haplotype_itr = std::find_if_not(std::next(haplotype_itr), std::cend(genotype), is_same_haplotype);
    }
    return true;
}

bool contains(const Genotype<Haplotype>& genotype, const Allele& allele);
bool includes(const Genotype<Haplotype>& genotype, const Allele& allele);

template <typename MappableType>
bool contains(const Genotype<MappableType>& genotype, const MappableType& element)
{
    return genotype.contains(element);
}

namespace detail {

template <typename MappableType>
bool contains(const Genotype<Haplotype>& lhs, const Genotype<MappableType>& rhs, std::true_type)
{
    using std::cbegin; using std::cend; using std::begin; using std::end;
    using AlleleReference = std::reference_wrapper<const MappableType>;
    if (lhs.ploidy() != rhs.ploidy()) return false;
    const auto lhs_copy = copy<MappableType>(lhs, mapped_region(rhs));
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
    using B = std::integral_constant<bool, std::is_same<MappableType1, Haplotype>::value
                                       && !std::is_same<MappableType2, Haplotype>::value>;
    return detail::contains(lhs, rhs, B {});
}

template <typename MappableType2, typename MappableType1>
bool are_equal_in_region(const Genotype<MappableType1>& lhs, const Genotype<MappableType1>& rhs,
                         const GenomicRegion& region)
{
    return copy<MappableType2>(lhs, region) == copy<MappableType2>(rhs, region);
}

template <typename MappableType>
bool is_homozygous(const Genotype<MappableType>& genotype, const MappableType& element)
{
    return genotype.count(element) == genotype.ploidy();
}

template <typename BinaryPredicate>
bool is_homozygous(const Genotype<Haplotype>& genotype, const Allele& allele, BinaryPredicate&& contains_pred)
{
    if (!contains(mapped_region(genotype), allele)) return true;
    if (genotype.is_homozygous()) return true;
    return all_of(genotype, [&] (const Haplotype& haplotype) { return contains_pred(haplotype, allele); });
}

template <typename BinaryPredicate>
bool is_heterozygous(const Genotype<Haplotype>& genotype, const Allele& allele, BinaryPredicate&& contains_pred)
{
    return !is_homozygous(genotype, allele, std::forward<BinaryPredicate>(contains_pred));
}

bool is_homozygous(const Genotype<Haplotype>& genotype, const Allele& allele);
bool is_heterozygous(const Genotype<Haplotype>& genotype, const Allele& allele);

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
        if (lhs.ploidy() < 2 || (lhs.ploidy() == 2 && !lhs.is_homozygous())) {
            return std::count_if(std::cbegin(lhs), std::cend(lhs),
                                 [&rhs] (const auto& element) { return rhs.contains(element); });
        } else {
            const auto& lhs_unique = lhs.copy_unique_ref();
            return std::count_if(std::cbegin(lhs_unique), std::cend(lhs_unique),
                                 [&rhs] (const auto& element) { return rhs.contains(element); });
        }
    } else {
        return count_shared(rhs, lhs);
    }
}

template <typename MappableType>
bool have_shared(const Genotype<MappableType>& lhs, const Genotype<MappableType>& rhs)
{
    if (lhs.ploidy() <= rhs.ploidy()) {
        return std::any_of(std::cbegin(lhs), std::cend(lhs), [&rhs] (const auto& element) { return rhs.contains(element); });
    } else {
        return have_shared(rhs, lhs);
    }
}

using GenotypeIndex = std::vector<unsigned>;

namespace detail {

namespace {

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
template <typename MappableType, typename Container>
auto construct_empty_genotype_container(const MappableBlock<MappableType, Container>& elements)
{
    return MappableBlock<Genotype<MappableType>, Container> {mapped_region(elements)};
}

} // namespace

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
auto generate_genotype(const Range& elements, const GenotypeIndex& element_indicies)
{
    GenotypeType<Range> result {static_cast<unsigned>(element_indicies.size())};
    for (const auto i : element_indicies) {
        result.emplace(elements[i]);
    }
    return result;
}

template <typename Range>
auto do_generate_all_genotypes(const Range& elements, const unsigned ploidy)
{
    using GenotypeTp = GenotypeType<Range>;
    using ResultType = decltype(construct_empty_genotype_container(elements));
    
    // Optimise simple cases
    if (ploidy == 0 || elements.empty()) {
        return construct_empty_genotype_container(elements);
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
    auto result = construct_empty_genotype_container(elements);
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

template <typename Range>
auto do_generate_all_genotypes(const Range& elements, const unsigned ploidy,
                               std::vector<GenotypeIndex>& indices)
{
    auto result = construct_empty_genotype_container(elements);
    if (ploidy == 0 || elements.empty()) {
        return result;
    }
    const auto num_elements = static_cast<unsigned>(elements.size());
    const auto result_size = num_genotypes(num_elements, ploidy);
    result.reserve(result_size);
    indices.reserve(result_size);
    GenotypeIndex element_indicies(ploidy, 0);
    while (true) {
        if (element_indicies[0] == num_elements) {
            unsigned i {0};
            while (++i < ploidy && element_indicies[i] == num_elements - 1);
            if (i == ploidy) break;
            ++element_indicies[i];
            std::fill_n(std::begin(element_indicies), i + 1, element_indicies[i]);
        }
        result.push_back(detail::generate_genotype(elements, element_indicies));
        indices.push_back(element_indicies);
        ++element_indicies[0];
    }
    return result;
}

template <typename Range, typename UnaryPredicate, typename OutputIterator>
auto do_generate_all_genotypes(const Range& elements, const unsigned ploidy,
                               UnaryPredicate pred, OutputIterator result_itr)
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

template <typename Range, typename UnaryPredicate, typename OutputIterator>
auto do_generate_all_genotypes(const Range& elements, const unsigned ploidy,
                               UnaryPredicate pred, OutputIterator result_itr,
                               std::vector<GenotypeIndex>& indices)
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
        if (pred(genotype)) {
            *result_itr++ = std::move(genotype);
            indices.push_back(element_indicies);
        }
        ++element_indicies[0];
    }
    return result_itr;
}

template <typename Range>
auto 
generate_all_genotypes(const Range& elements, const unsigned ploidy,
                       std::true_type)
{
    using MappableType = value_type_t<typename Range::value_type>;
    std::vector<std::shared_ptr<MappableType>> temp_pointers(elements.size());
    std::transform(std::cbegin(elements), std::cend(elements), std::begin(temp_pointers),
                   [] (const MappableType& element) { return std::make_shared<MappableType>(element); });
    return do_generate_all_genotypes(temp_pointers, ploidy);
}

template <typename Range>
auto 
generate_all_genotypes(const Range& elements, const unsigned ploidy,
                       std::vector<GenotypeIndex>& indices, std::true_type)
{
    using MappableType = value_type_t<typename Range::value_type>;
    std::vector<std::shared_ptr<MappableType>> temp_pointers(elements.size());
    std::transform(std::cbegin(elements), std::cend(elements), std::begin(temp_pointers),
                   [] (const MappableType& element) { return std::make_shared<MappableType>(element); });
    return do_generate_all_genotypes(temp_pointers, ploidy, indices);
}

template <typename Range>
auto 
generate_all_genotypes(const Range& elements, const unsigned ploidy,
                       std::false_type)
{
    return do_generate_all_genotypes(elements, ploidy);
}

template <typename Range>
auto 
generate_all_genotypes(const Range& elements, const unsigned ploidy,
                       std::vector<GenotypeIndex>& indices, std::false_type)
{
    return do_generate_all_genotypes(elements, ploidy, indices);
}

template <typename Range, typename UnaryPredicate, typename OutputIterator>
OutputIterator
generate_all_genotypes(const Range& elements, const unsigned ploidy,
                       UnaryPredicate selector, OutputIterator result_itr,
                       std::true_type)
{
    using MappableType = value_type_t<typename Range::value_type>;
    std::vector<std::shared_ptr<MappableType>> temp_pointers(elements.size());
    std::transform(std::cbegin(elements), std::cend(elements), std::begin(temp_pointers),
                   [] (const MappableType& element) { return std::make_shared<MappableType>(element); });
    return do_generate_all_genotypes(temp_pointers, ploidy, selector, result_itr);
}

template <typename Range, typename UnaryPredicate, typename OutputIterator>
OutputIterator
generate_all_genotypes(const Range& elements, const unsigned ploidy,
                       UnaryPredicate selector, OutputIterator result_itr, std::vector<GenotypeIndex>& indices,
                       std::true_type)
{
    using MappableType = value_type_t<typename Range::value_type>;
    std::vector<std::shared_ptr<MappableType>> temp_pointers(elements.size());
    std::transform(std::cbegin(elements), std::cend(elements), std::begin(temp_pointers),
                   [] (const MappableType& element) { return std::make_shared<MappableType>(element); });
    return do_generate_all_genotypes(temp_pointers, ploidy, selector, result_itr, indices);
}

template <typename Range, typename UnaryPredicate, typename OutputIterator>
OutputIterator
generate_all_genotypes(const Range& elements, const unsigned ploidy,
                       UnaryPredicate selector, OutputIterator result_itr,
                       std::false_type)
{
    return do_generate_all_genotypes(elements, ploidy, selector, result_itr);
}

template <typename Range, typename UnaryPredicate, typename OutputIterator>
OutputIterator
generate_all_genotypes(const Range& elements, const unsigned ploidy,
                       UnaryPredicate selector, OutputIterator result_itr, std::vector<GenotypeIndex>& indices,
                       std::false_type)
{
    return do_generate_all_genotypes(elements, ploidy, selector, result_itr, indices);
}

} // namespace detail

template <typename Range>
auto
generate_all_genotypes(const Range& elements, const unsigned ploidy)
{
    using MappableType = detail::value_type_t<typename Range::value_type>;
    return detail::generate_all_genotypes(elements, ploidy, detail::share_memory<MappableType> {});
}

template <typename Range>
auto
generate_all_genotypes(const Range& elements, const unsigned ploidy,
                       std::vector<GenotypeIndex>& indices)
{
    using MappableType = detail::value_type_t<typename Range::value_type>;
    return detail::generate_all_genotypes(elements, ploidy, indices, detail::share_memory<MappableType> {});
}

template <typename Range, typename UnaryPredicate, typename OutputIterator>
OutputIterator
generate_all_genotypes(const Range& elements, const unsigned ploidy,
                       UnaryPredicate selector, OutputIterator result_itr)
{
    using MappableType = detail::value_type_t<typename Range::value_type>;
    return detail::generate_all_genotypes(elements, ploidy, selector, result_itr, detail::share_memory<MappableType> {});
}

template <typename Range, typename UnaryPredicate, typename OutputIterator>
OutputIterator
generate_all_genotypes(const Range& elements, const unsigned ploidy,
                       UnaryPredicate selector, OutputIterator result_itr, std::vector<GenotypeIndex>& indices)
{
    using MappableType = detail::value_type_t<typename Range::value_type>;
    return detail::generate_all_genotypes(elements, ploidy, selector, result_itr, indices, detail::share_memory<MappableType> {});
}

std::vector<Genotype<Haplotype>>
generate_all_genotypes(const std::vector<std::shared_ptr<Haplotype>>& haplotypes, unsigned ploidy);

template <typename MappableType>
bool is_max_zygosity(const Genotype<MappableType>& genotype)
{
    return genotype.zygosity() == genotype.ploidy();
}

std::size_t num_max_zygosity_genotypes(unsigned num_elements, unsigned ploidy);
boost::optional<std::size_t> num_max_zygosity_genotypes_noexcept(unsigned num_elements, unsigned ploidy) noexcept;

template <typename Range>
void
generate_all_max_zygosity_genotypes(std::vector<Genotype<detail::value_type_t<typename Range::value_type>>>& result,
                                    const Range& elements, const unsigned ploidy)
{
    if (elements.size() >= ploidy) {
        try {
            result.reserve(num_max_zygosity_genotypes(elements.size(), ploidy));
        } catch (const std::overflow_error& e) {
            result.reserve(std::numeric_limits<std::size_t>::max()); // Probably this will throw a bad_alloc
        }
        generate_all_genotypes(elements, ploidy, [ploidy] (const auto& genotype) { return genotype.zygosity() == ploidy; },
                               std::back_inserter(result));
    }
}

template <typename Range>
auto
generate_all_max_zygosity_genotypes(const Range& elements, const unsigned ploidy)
{
    using MappableType = detail::value_type_t<typename Range::value_type>;
    std::vector<Genotype<MappableType>> result {};
    generate_all_max_zygosity_genotypes(result, elements, ploidy);
    return result;
}

template <typename Range>
void
generate_all_max_zygosity_genotypes(std::vector<Genotype<detail::value_type_t<typename Range::value_type>>>& result,
                                    std::vector<GenotypeIndex>& indices,
                                    const Range& elements, const unsigned ploidy)
{
    if (elements.size() >= ploidy) {
        try {
            result.reserve(num_max_zygosity_genotypes(elements.size(), ploidy));
        } catch (const std::overflow_error& e) {
            result.reserve(std::numeric_limits<std::size_t>::max()); // Probably this will throw a bad_alloc
        }
        generate_all_genotypes(elements, ploidy, [ploidy] (const auto& genotype) { return genotype.zygosity() == ploidy; },
                               std::back_inserter(result), indices);
    }
}

template <typename Range>
auto
generate_all_max_zygosity_genotypes(const Range& elements, const unsigned ploidy,
                                    std::vector<GenotypeIndex>& indices)
{
    using MappableType = detail::value_type_t<typename Range::value_type>;
    std::vector<Genotype<MappableType>> result {};
    generate_all_max_zygosity_genotypes(result, indices, elements, ploidy);
    return result;
}

template <typename UnaryPredicate>
std::vector<Genotype<Haplotype>>
extend_genotypes(const std::vector<Genotype<Haplotype>>& genotypes,
                 const std::vector<Haplotype>& haplotypes,
                 UnaryPredicate selector)
{
    std::vector<std::shared_ptr<Haplotype>> temp_pointers(haplotypes.size());
    std::transform(std::cbegin(haplotypes), std::cend(haplotypes), std::begin(temp_pointers),
                   [] (const auto& haplotype) { return std::make_shared<Haplotype>(haplotype); });
    std::vector<Genotype<Haplotype>> result {};
    result.reserve(genotypes.size() * haplotypes.size());
    for (const auto& genotype : genotypes) {
        for (const auto& haplotype_ptr : temp_pointers) {
            if (selector(genotype, *haplotype_ptr)) {
                auto extended_genotype = genotype;
                extended_genotype.emplace(haplotype_ptr);
                result.push_back(std::move(extended_genotype));
            }
        }
    }
    return result;
}

template <typename UnaryPredicate>
std::pair<std::vector<Genotype<Haplotype>>, std::vector<GenotypeIndex>>
extend_genotypes(const std::vector<Genotype<Haplotype>>& genotypes,
                 const std::vector<GenotypeIndex>& indices,
                 const std::vector<Haplotype>& haplotypes,
                 UnaryPredicate selector)
{
    std::vector<std::shared_ptr<Haplotype>> temp_pointers(haplotypes.size());
    std::transform(std::cbegin(haplotypes), std::cend(haplotypes), std::begin(temp_pointers),
                   [] (const auto& haplotype) { return std::make_shared<Haplotype>(haplotype); });
    std::vector<Genotype<Haplotype>> extended_genotypes {};
    std::vector<GenotypeIndex> extended_indices {};
    const auto max_extended_genotypes = genotypes.size() * haplotypes.size();
    extended_genotypes.reserve(max_extended_genotypes);
    extended_indices.reserve(max_extended_genotypes);
    for (std::size_t g {0}; g < genotypes.size(); ++g) {
        for (std::size_t h {0}; h < haplotypes.size(); ++h) {
            const auto& haplotype_ptr = temp_pointers[h];
            if (selector(genotypes[g], haplotypes[h])) {
                auto extended_genotype = genotypes[g];
                extended_genotype.emplace(haplotype_ptr);
                extended_genotypes.push_back(std::move(extended_genotype));
                auto new_index = indices[g];
                new_index.push_back(h);
                extended_indices.push_back(std::move(new_index));
            }
        }
    }
    return std::make_pair(std::move(extended_genotypes), std::move(extended_indices));
}

inline
std::pair<std::vector<Genotype<Haplotype>>, std::vector<GenotypeIndex>>
extend_genotypes(const std::vector<Genotype<Haplotype>>& genotypes,
                 const std::vector<GenotypeIndex>& indices,
                 const std::vector<Haplotype>& haplotypes)
{
    const static auto default_selector = [] (const auto&, const auto&) noexcept { return true; };
    return extend_genotypes(genotypes, indices, haplotypes, default_selector);
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
    result.reserve(genotype.zygosity());
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

template <typename S>
void print_alleles(S&& stream, const Genotype<Haplotype>& genotype)
{
    if (genotype.ploidy() == 0) {
        stream << "[]";
    }
    const auto haplotype_counts = make_element_count_map(genotype);
    std::vector<std::pair<Haplotype, unsigned>> p {haplotype_counts.begin(), haplotype_counts.end()};
    stream << "[";
    for (unsigned i {0}; i < p.size() - 1; ++i) {
        print_alleles(stream, p[i].first);
        stream << "(" << p[i].second << "),";
    }
    print_alleles(stream, p.back().first);
    stream << "(" << p.back().second << ")]";
}

template <typename S>
void print_variant_alleles(S&& stream, const Genotype<Haplotype>& genotype)
{
    if (genotype.ploidy() == 0) {
        stream << "[]";
    }
    
    const auto unique_haplotypes = genotype.copy_unique();
    
    stream << "[";
    for (unsigned i {0}; i < unique_haplotypes.size() - 1; ++i) {
        print_variant_alleles(stream, unique_haplotypes[i]);
        stream << "(" << genotype.count(unique_haplotypes[i]) << "),";
    }
    
    print_variant_alleles(stream, unique_haplotypes.back());
    stream << "(" << genotype.count(unique_haplotypes.back()) << ")]";
}

void print_alleles(const Genotype<Haplotype>& genotype);
void print_variant_alleles(const Genotype<Haplotype>& genotype);

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
