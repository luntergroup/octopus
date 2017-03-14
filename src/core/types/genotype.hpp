// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef genotype_hpp
#define genotype_hpp

#include <vector>
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

#include "concepts/equitable.hpp"
#include "concepts/mappable.hpp"
#include "allele.hpp"
#include "haplotype.hpp"

namespace octopus {

template <typename T>
using EnableIfGenotypable = std::enable_if_t<
                                std::is_same<T, Haplotype>::value
                                || std::is_same<T, Allele>::value
                                || std::is_same<T, ContigAllele>::value
                            >;

template <typename MappableType, typename = EnableIfGenotypable<MappableType>> class Genotype;

template <typename MappableType>
class Genotype<MappableType> : public Equitable<Genotype<MappableType>>, public Mappable<Genotype<MappableType>>
{
public:
    using ElementType   = MappableType;
    using MappingDomain = RegionType<ElementType>;
    
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
    
    std::vector<MappableType> copy_unique() const;
    
private:
    std::vector<MappableType> elements_;

public:
    using Iterator = typename decltype(elements_)::const_iterator;
    
    Iterator begin() const noexcept;
    Iterator end() const noexcept ;
    Iterator cbegin() const noexcept ;
    Iterator cend() const noexcept ;
};

template <>
class Genotype<Haplotype> : public Equitable<Genotype<Haplotype>>, public Mappable<Genotype<Haplotype>>
{
public:
    using ElementType   = Haplotype;
    using MappingDomain = Haplotype::MappingDomain;
    
    Genotype() = default;
    
    explicit Genotype(unsigned ploidy);
    explicit Genotype(unsigned ploidy, const Haplotype& init);
    explicit Genotype(unsigned ploidy, const std::shared_ptr<Haplotype>& init);
    explicit Genotype(std::initializer_list<Haplotype> elements);
    explicit Genotype(std::initializer_list<std::shared_ptr<Haplotype>> elements);
    
    Genotype(const Genotype&)            = default;
    Genotype& operator=(const Genotype&) = default;
    Genotype(Genotype&&)                 = default;
    Genotype& operator=(Genotype&&)      = default;
    
    ~Genotype() = default;
    
    template <typename T> void emplace(T&& element);
    void emplace(const std::shared_ptr<Haplotype>& element);
    
    const Haplotype& operator[](unsigned n) const;
    
    const GenomicRegion& mapped_region() const noexcept;
    
    unsigned ploidy() const noexcept;
    
    bool contains(const Haplotype& haplotype) const;
    unsigned count(const Haplotype& haplotype) const;
    
    bool is_homozygous() const;
    unsigned zygosity() const;
    
    std::vector<Haplotype> copy_unique() const;
    std::vector<std::reference_wrapper<const Haplotype>> copy_unique_ref() const;
    
private:
    using HaplotypePtr  = std::shared_ptr<Haplotype>;
    using BaseContainer = std::vector<HaplotypePtr>;
    using BaseIterator  = typename BaseContainer::const_iterator;
    
    BaseContainer haplotypes_;
    
    struct HaplotypePtrLess
    {
        bool operator()(const HaplotypePtr& lhs, const HaplotypePtr& rhs) const;
        bool operator()(const Haplotype& lhs, const HaplotypePtr& rhs) const;
        bool operator()(const HaplotypePtr& lhs, const Haplotype& rhs) const;
    };
    
    struct HaplotypePtrEqual
    {
        bool operator()(const HaplotypePtr& lhs, const HaplotypePtr& rhs) const;
        bool operator()(const Haplotype& lhs, const HaplotypePtr& rhs) const;
        bool operator()(const HaplotypePtr& lhs, const Haplotype& rhs) const;
    };
    
public:
    class Iterator : public BaseIterator
    {
    public:
        using value_type = Haplotype;
        using reference  = const Haplotype&;
        using pointer    = const Haplotype*;
        
        Iterator(BaseIterator it);
        reference operator*() const;
    };
    
    Iterator begin() const noexcept;
    Iterator end() const noexcept;
    Iterator cbegin() const noexcept;
    Iterator cend() const noexcept;
};

// Genotype<MappableType>

template <typename MappableType>
template <typename T>
void Genotype<MappableType>::emplace(T&& element)
{
    elements_.emplace_back(std::forward<T>(element));
}

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
{}

template <typename MappableType>
const typename Genotype<MappableType>::MappingDomain& Genotype<MappableType>::mapped_region() const noexcept
{
    return elements_.front().mapped_region();
}

template <typename MappableType>
const MappableType& Genotype<MappableType>::operator[](const unsigned n) const
{
    return elements_[n];
}

template <typename MappableType>
unsigned Genotype<MappableType>::ploidy() const noexcept
{
    return static_cast<unsigned>(elements_.size());
}

template <typename MappableType>
bool Genotype<MappableType>::is_homozygous() const
{
    return std::adjacent_find(std::cbegin(elements_), std::cend(elements_),
                              std::not_equal_to<Allele>()) == std::cend(elements_);
}

template <typename MappableType>
unsigned Genotype<MappableType>::zygosity() const
{
    if (ploidy() == 1 || is_homozygous()) {
        return 1;
    } else if (ploidy() == 2) {
        return 2;
    }
    return static_cast<unsigned>(copy_unique().size());
}

template <typename MappableType>
bool Genotype<MappableType>::contains(const MappableType& element) const
{
    return std::find(std::cbegin(elements_), std::cend(elements_), element) != std::cend(elements_);
}

template <typename MappableType>
unsigned Genotype<MappableType>::count(const MappableType& element) const
{
    return static_cast<unsigned>(std::count(std::cbegin(elements_), std::cend(elements_), element));
}

template <typename MappableType>
std::vector<MappableType> Genotype<MappableType>::copy_unique() const
{
    auto result = elements_;
    std::sort(std::begin(result), std::end(result));
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

template <typename MappableType>
typename Genotype<MappableType>::Iterator Genotype<MappableType>::begin() const noexcept
{
    return std::begin(elements_);
}

template <typename MappableType>
typename Genotype<MappableType>::Iterator Genotype<MappableType>::end() const noexcept
{
    return std::end(elements_);
}

template <typename MappableType>
typename Genotype<MappableType>::Iterator Genotype<MappableType>::cbegin() const noexcept
{
    return std::cbegin(elements_);
}

template <typename MappableType>
typename Genotype<MappableType>::Iterator Genotype<MappableType>::cend() const noexcept
{
    return std::cend(elements_);
}

// Genotype<Haplotype>

template <typename T>
void Genotype<Haplotype>::emplace(T&& haplotype)
{
    haplotypes_.emplace_back(std::make_shared<Haplotype>(std::forward<T>(haplotype)));
    std::inplace_merge(std::begin(haplotypes_), std::prev(std::end(haplotypes_)),
                       std::end(haplotypes_), HaplotypePtrLess {});
}

// non-member methods

template <typename MappableType>
bool is_haploid(const Genotype<MappableType>& genotype)
{
    return genotype.ploidy() == 1;
}

template <typename MappableType>
bool is_diploid(const Genotype<MappableType>& genotype)
{
    return genotype.ploidy() == 2;
}

template <typename MappableType>
bool is_triploid(const Genotype<MappableType>& genotype)
{
    return genotype.ploidy() == 3;
}

template <typename MappableType>
bool is_tetraploid(const Genotype<MappableType>& genotype)
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

template <typename MappableType>
struct ShareMemory : public std::is_same<MappableType, Haplotype> {};

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
};

template <typename MappableType, typename Container>
auto copy_each_cached(const Container& genotypes, const GenomicRegion& region)
{
    using MappableType2 = copy_type<typename Container::value_type>;
    return copy_each_cached<MappableType>(genotypes, region, ShareMemory<MappableType2> {});
}

} // namespace detail

template <typename MappableType, typename Container>
std::vector<Genotype<MappableType>> copy_each(const Container& genotypes, const GenomicRegion& region)
{
    if (genotypes.size() < 10) {
        return detail::copy_each_basic<MappableType>(genotypes, region);
    } else {
        return detail::copy_each_cached<MappableType>(genotypes, region);
    }
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

bool is_homozygous(const Genotype<Haplotype>& genotype, const Allele& allele);

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
std::size_t element_cardinality_in_genotypes(unsigned num_elements, unsigned ploidy);

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

template <typename Container>
using GenotypeType = Genotype<typename ValueType<typename Container::value_type>::type>;

} // namespace

template <typename Container>
auto generate_all_haploid_genotypes(const Container& elements)
{
    using ResultType = std::vector<GenotypeType<Container>>;
    ResultType result {};
    result.reserve(elements.size());
    for (const auto& element : elements) {
        result.emplace_back(1, element);
    }
    return result;
}

template <typename Container>
auto generate_all_diploid_biallelic_genotypes(const Container& elements)
{
    using GenotypeTp = GenotypeType<Container>;
    using ResultType = std::vector<GenotypeTp>;
    return ResultType {
        GenotypeTp {elements[0], elements[0]},
        GenotypeTp {elements[0], elements[1]},
        GenotypeTp {elements[1], elements[1]}
    };
}

template <typename Container>
auto generate_all_diploid_triallelic_genotypes(const Container& elements)
{
    using GenotypeTp = GenotypeType<Container>;
    using ResultType = std::vector<GenotypeTp>;
    return ResultType {
        GenotypeTp {elements[0], elements[0]},
        GenotypeTp {elements[0], elements[1]},
        GenotypeTp {elements[0], elements[2]},
        GenotypeTp {elements[1], elements[1]},
        GenotypeTp {elements[1], elements[2]},
        GenotypeTp {elements[2], elements[2]}
    };
}

template <typename Container>
auto generate_all_diploid_tetraallelic_genotypes(const Container& elements)
{
    using GenotypeTp = GenotypeType<Container>;
    using ResultType = std::vector<GenotypeTp>;
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

template <typename Container>
auto generate_all_triploid_biallelic_genotypes(const Container& elements)
{
    using GenotypeTp = GenotypeType<Container>;
    using ResultType = std::vector<GenotypeTp>;
    return ResultType {
        GenotypeTp {elements[0], elements[0], elements[0]},
        GenotypeTp {elements[0], elements[0], elements[1]},
        GenotypeTp {elements[0], elements[1], elements[1]},
        GenotypeTp {elements[1], elements[1], elements[1]}
    };
}

template <typename Container>
auto generate_genotype(const Container& elements, const std::vector<unsigned>& element_indicies)
{
    GenotypeType<Container> result {static_cast<unsigned>(element_indicies.size())};
    for (const auto i : element_indicies) {
        result.emplace(elements[i]);
    }
    return result;
}

template <typename Container>
auto do_generate_all_genotypes(const Container& elements, const unsigned ploidy)
{
    using GenotypeTp = GenotypeType<Container>;
    using ResultType = std::vector<GenotypeTp>;
    
    // Optimise simple cases
    if (ploidy == 0 || elements.empty()) {
        return ResultType {};
    }
    const auto num_elements = static_cast<unsigned>(elements.size());
    if (num_elements == 1) {
        return ResultType {GenotypeTp {ploidy, elements.front()}};
    }
    if (ploidy == 2) {
        switch(num_elements) {
            case 2: return detail::generate_all_diploid_biallelic_genotypes(elements);
            case 3: return detail::generate_all_diploid_triallelic_genotypes(elements);
            case 4: return detail::generate_all_diploid_tetraallelic_genotypes(elements);
        }
    }
    if (ploidy == 1) {
        return detail::generate_all_haploid_genotypes(elements);
    }
    if (ploidy == 3 && num_elements == 2) {
        return detail::generate_all_triploid_biallelic_genotypes(elements);
    }
    
    // Otherwise resort to general algorithm
    ResultType result {};
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

template <typename MappableType>
struct RequiresSharedMemory : public std::is_same<MappableType, Haplotype> {};

template <typename MappableType>
auto generate_all_genotypes(const std::vector<MappableType>& elements, const unsigned ploidy,
                            std::true_type)
{
    std::vector<std::shared_ptr<MappableType>> temp_pointers(elements.size());
    std::transform(std::cbegin(elements), std::cend(elements), std::begin(temp_pointers),
                   [] (const auto& element) {
                       return std::make_shared<MappableType>(element);
                   });
    return do_generate_all_genotypes(temp_pointers, ploidy);
}

template <typename MappableType>
auto generate_all_genotypes(const std::vector<std::reference_wrapper<const MappableType>>& elements,
                            const unsigned ploidy,
                            std::true_type)
{
    std::vector<std::shared_ptr<MappableType>> temp_pointers(elements.size());
    std::transform(std::cbegin(elements), std::cend(elements), std::begin(temp_pointers),
                   [] (const auto& element) {
                       return std::make_shared<MappableType>(element.get());
                   });
    return do_generate_all_genotypes(temp_pointers, ploidy);
}

template <typename MappableType>
auto generate_all_genotypes(const std::vector<MappableType>& elements, const unsigned ploidy,
                            std::false_type)
{
    return do_generate_all_genotypes(elements, ploidy);
}

} // namespace detail

template <typename MappableType>
std::vector<Genotype<MappableType>>
generate_all_genotypes(const std::vector<MappableType>& elements, const unsigned ploidy)
{
    return detail::generate_all_genotypes(elements, ploidy, detail::RequiresSharedMemory<MappableType> {});
}

template <typename MappableType>
std::vector<Genotype<MappableType>>
generate_all_genotypes(const std::vector<std::reference_wrapper<const MappableType>>& elements,
                       const unsigned ploidy)
{
    return detail::generate_all_genotypes(elements, ploidy, detail::RequiresSharedMemory<MappableType> {});
}

std::vector<Genotype<Haplotype>>
generate_all_genotypes(const std::vector<std::shared_ptr<Haplotype>>& haplotypes, unsigned ploidy);

namespace detail {

inline std::size_t estimate_num_elements(const std::size_t num_genotypes)
{
    return num_genotypes;
}

} // namespace detail

template <typename MappableType>
auto extract_all_elements(const std::vector<Genotype<MappableType>>& genotypes)
{
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

template <typename MappableType>
auto extract_all_element_refs(const std::vector<Genotype<MappableType>>& genotypes)
{
    std::unordered_set<std::reference_wrapper<const MappableType>> unique_elements {};
    unique_elements.reserve();
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
