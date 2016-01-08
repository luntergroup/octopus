//
//  genotype.hpp
//  Octopus
//
//  Created by Daniel Cooke on 24/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__genotype__
#define __Octopus__genotype__

#include <vector>
#include <unordered_map>
#include <initializer_list>
#include <ostream>
#include <cstddef>
#include <iterator>
#include <algorithm>
#include <type_traits>

#include <boost/functional/hash.hpp>

#include "allele.hpp"
#include "haplotype.hpp"
#include "equitable.hpp"
#include "mappable.hpp"

// Genotype should only store Haplotype and Alleles
template <typename MappableType, typename = std::enable_if_t<std::is_base_of<Mappable<MappableType>, MappableType>::value>>
class Genotype;

template <typename MappableType>
class Genotype<MappableType> : public Equitable<Genotype<MappableType>>, public Mappable<Genotype<MappableType>>
{
public:
    using Iterator = typename std::vector<MappableType>::const_iterator;
    
    Genotype() = default;
    explicit Genotype(unsigned ploidy);
    explicit Genotype(unsigned ploidy, const MappableType& init);
    explicit Genotype(std::initializer_list<MappableType> elements);
    ~Genotype() = default;
    
    Genotype(const Genotype&)            = default;
    Genotype& operator=(const Genotype&) = default;
    Genotype(Genotype&&)                 = default;
    Genotype& operator=(Genotype&&)      = default;
    
    const MappableType& at(unsigned n) const;
    const MappableType& operator[](unsigned n) const;
    template <typename T> void emplace(T&& element);
    
    Iterator begin() const noexcept;
    Iterator end() const noexcept ;
    Iterator cbegin() const noexcept ;
    Iterator cend() const noexcept ;
    
    const GenomicRegion& get_region() const noexcept;
    
    unsigned ploidy() const noexcept;
    
    bool contains(const MappableType& element) const;
    unsigned count(const MappableType& element) const;
    bool is_homozygous() const;
    unsigned zygosity() const;
    
    std::vector<MappableType> copy_unique() const;
    
private:
    std::vector<MappableType> elements_;
};

template <typename MappableType>
Genotype<MappableType>::Genotype(const unsigned ploidy)
:
elements_ {}
{
    elements_.reserve(ploidy);
}

template <typename MappableType>
Genotype<MappableType>::Genotype(const unsigned ploidy, const MappableType& init)
:
elements_ {ploidy, init}
{}

template <typename MappableType>
Genotype<MappableType>::Genotype(std::initializer_list<MappableType> elements)
:
elements_ {elements}
{
    std::sort(elements_.begin(), elements_.end());
}

template <typename MappableType>
const MappableType& Genotype<MappableType>::at(const unsigned n) const
{
    return elements_.at(n);
}

template <typename MappableType>
const MappableType& Genotype<MappableType>::operator[](const unsigned n) const
{
    return elements_[n];
}

template <typename MappableType>
template <typename T>
void Genotype<MappableType>::emplace(T&& element)
{
    elements_.emplace_back(std::forward<T>(element));
    std::inplace_merge(std::begin(elements_), std::prev(std::end(elements_)), std::end(elements_));
}

template <typename MappableType>
typename Genotype<MappableType>::Iterator Genotype<MappableType>::begin() const noexcept
{
    return elements_.begin();
}

template <typename MappableType>
typename Genotype<MappableType>::Iterator Genotype<MappableType>::end() const noexcept
{
    return elements_.end();
}

template <typename MappableType>
typename Genotype<MappableType>::Iterator Genotype<MappableType>::cbegin() const noexcept
{
    return elements_.cbegin();
}

template <typename MappableType>
typename Genotype<MappableType>::Iterator Genotype<MappableType>::cend() const noexcept
{
    return elements_.cend();
}

template <typename MappableType>
const GenomicRegion& Genotype<MappableType>::get_region() const noexcept
{
    return elements_.front().get_region();
}

template <typename MappableType>
unsigned Genotype<MappableType>::ploidy() const noexcept
{
    return static_cast<unsigned>(elements_.size());
}

template <typename MappableType>
bool Genotype<MappableType>::is_homozygous() const
{
    return elements_.front() == elements_.back();
}

template <typename MappableType>
unsigned Genotype<MappableType>::zygosity() const
{
    unsigned result {};
    
    for (auto first = std::cbegin(elements_), last = std::cend(elements_); first != last;) {
        ++result;
        first = std::upper_bound(first, last, *first);
    }
    
    return result;
}

template <typename MappableType>
bool Genotype<MappableType>::contains(const MappableType& element) const
{
    return std::binary_search(std::cbegin(elements_), std::cend(elements_), element);
}

template <typename MappableType>
unsigned Genotype<MappableType>::count(const MappableType& element) const
{
    const auto equal_range = std::equal_range(std::cbegin(elements_), std::cend(elements_), element);
    return static_cast<unsigned>(std::distance(equal_range.first, equal_range.second));
}

template <typename MappableType>
std::vector<MappableType> Genotype<MappableType>::copy_unique() const
{
    std::vector<MappableType> result {};
    result.reserve(ploidy());
    
    std::unique_copy(std::cbegin(elements_), std::cend(elements_), std::back_inserter(result));
    
    result.shrink_to_fit();
    
    return result;
}

template <>
class Genotype<Allele> : public Equitable<Genotype<Allele>>, public Mappable<Genotype<Allele>>
{
public:
    using Iterator = typename std::vector<Allele>::const_iterator;
    
    Genotype() = default;
    explicit Genotype(unsigned ploidy);
    explicit Genotype(unsigned ploidy, const Allele& init);
    explicit Genotype(std::initializer_list<Allele> alleles);
    ~Genotype() = default;
    
    Genotype(const Genotype&)            = default;
    Genotype& operator=(const Genotype&) = default;
    Genotype(Genotype&&)                 = default;
    Genotype& operator=(Genotype&&)      = default;
    
    const Allele& at(unsigned n) const;
    const Allele& operator[](unsigned n) const;
    template <typename T> void emplace(T&& allele);
    
    Iterator begin() const noexcept;
    Iterator end() const noexcept ;
    Iterator cbegin() const noexcept ;
    Iterator cend() const noexcept ;
    
    unsigned ploidy() const noexcept;
    
    const GenomicRegion& get_region() const noexcept;
    
    bool contains(const Allele& allele) const;
    unsigned count(const Allele& allele) const;
    bool is_homozygous() const;
    unsigned zygosity() const;
    
    std::vector<Allele> copy_unique() const;
    
private:
    std::vector<Allele> alleles_;
};

template <typename T>
void Genotype<Allele>::emplace(T&& allele)
{
    alleles_.emplace_back(std::forward<T>(allele));
}

// non-member methods

template <typename MappableType2, typename MappableType1>
Genotype<MappableType2> splice(const Genotype<MappableType1>& genotype, const GenomicRegion& region)
{
    Genotype<MappableType2> result {genotype.ploidy()};
    
    for (const auto& mappable : genotype) {
        result.emplace(splice<MappableType2>(mappable, region));
    }
    
    return result;
}

bool contains(const Genotype<Haplotype>& genotype, const Allele& allele);
bool contains_exact(const Genotype<Haplotype>& genotype, const Allele& allele);

template <typename MappableType>
bool contains(const Genotype<MappableType>& genotype, const MappableType& element)
{
    return genotype.contains(element);
}

template <typename MappableType1, typename MappableType2>
bool contains(const Genotype<MappableType1>& lhs, const Genotype<MappableType2>& rhs)
{
    return splice<MappableType2>(lhs, rhs[0].get_region()) == rhs;
}

template <typename MappableType>
bool are_equal_in_region(const Genotype<MappableType>& lhs, const Genotype<MappableType>& rhs,
                         const GenomicRegion& region)
{
    return splice<MappableType>(lhs, region) == splice<MappableType>(rhs, region);
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

namespace std
{
    template <typename MappableType> struct hash<Genotype<MappableType>>
    {
        size_t operator()(const Genotype<MappableType>& genotype) const
        {
            return boost::hash_range(std::cbegin(genotype), std::cend(genotype));
        }
    };
    
    template <typename MappableType> struct hash<reference_wrapper<const Genotype<MappableType>>>
    {
        size_t operator()(const reference_wrapper<const Genotype<MappableType>> genotype) const
        {
            return hash<Genotype<MappableType>>()(genotype);
        }
    };
} // namespace std

size_t num_genotypes(unsigned num_elements, unsigned ploidy);

namespace detail
{
    template <typename MappableType>
    Genotype<MappableType> generate_genotype(const std::vector<MappableType>& elements,
                                             const std::vector<unsigned>& element_indicies)
    {
        Genotype<MappableType> result {static_cast<unsigned>(element_indicies.size())};
        
        for (auto i : element_indicies) {
            result.emplace(elements[i]);
        }
        
        return result;
    }
    
    template <typename MappableType>
    std::vector<Genotype<MappableType>>
    generate_all_haploid_genotypes(const std::vector<MappableType>& elements)
    {
        std::vector<Genotype<MappableType>> result {};
        result.reserve(elements.size());
        
        for (const auto& element : elements) {
            result.emplace_back(1, element);
        }
        
        return result;
    }
    
    template <typename MappableType>
    std::vector<Genotype<MappableType>>
    generate_all_diploid_biallelic_genotypes(const std::vector<MappableType>& elements)
    {
        return std::vector<Genotype<MappableType>> {
            Genotype<MappableType> {elements.front(), elements.front()},
            Genotype<MappableType> {elements.front(), elements.back()},
            Genotype<MappableType> {elements.back(), elements.back()}
        };
    }
    
    template <typename MappableType>
    std::vector<Genotype<MappableType>>
    generate_all_diploid_triallelic_genotypes(const std::vector<MappableType>& elements)
    {
        return std::vector<Genotype<MappableType>> {
            Genotype<MappableType> {elements[0], elements[0]},
            Genotype<MappableType> {elements[0], elements[1]},
            Genotype<MappableType> {elements[0], elements[2]},
            Genotype<MappableType> {elements[1], elements[1]},
            Genotype<MappableType> {elements[1], elements[2]},
            Genotype<MappableType> {elements[2], elements[2]}
        };
    }
    
    template <typename MappableType>
    std::vector<Genotype<MappableType>>
    generate_all_diploid_tetraallelic_genotypes(const std::vector<MappableType>& elements)
    {
        return std::vector<Genotype<MappableType>> {
            Genotype<MappableType> {elements[0], elements[0]},
            Genotype<MappableType> {elements[0], elements[1]},
            Genotype<MappableType> {elements[0], elements[2]},
            Genotype<MappableType> {elements[0], elements[3]},
            Genotype<MappableType> {elements[1], elements[1]},
            Genotype<MappableType> {elements[1], elements[2]},
            Genotype<MappableType> {elements[1], elements[3]},
            Genotype<MappableType> {elements[2], elements[2]},
            Genotype<MappableType> {elements[2], elements[3]},
            Genotype<MappableType> {elements[3], elements[3]}
        };
    }
} // namespace detail

// Assumes the input Haplotypes are unique
template <typename MappableType>
std::vector<Genotype<MappableType>>
generate_all_genotypes(const std::vector<MappableType>& elements, const unsigned ploidy)
{
    // simple cases are optimised
    
    if (ploidy == 0 || elements.empty()) return {};
    
    const auto num_elements = static_cast<unsigned>(elements.size());
    
    if (num_elements == 1) return {Genotype<MappableType> {ploidy, elements.front()}};
    
    if (ploidy == 2) {
        if (num_elements == 2) return detail::generate_all_diploid_biallelic_genotypes(elements);
        if (num_elements == 3) return detail::generate_all_diploid_triallelic_genotypes(elements);
        if (num_elements == 4) return detail::generate_all_diploid_tetraallelic_genotypes(elements);
    }
    
    if (ploidy == 1) return detail::generate_all_haploid_genotypes(elements);
    
    std::vector<Genotype<MappableType>> result {};
    result.reserve(num_genotypes(num_elements, ploidy));
    
    std::vector<unsigned> element_indicies(ploidy, 0);
    
    unsigned i {};
    
    while (true) {
        if (element_indicies[i] == num_elements) {
            while (i++ < ploidy && element_indicies[i] == num_elements - 1);
            
            if (i == ploidy) break;
            
            ++element_indicies[i];
            
            std::fill_n(std::begin(element_indicies), i + 1, element_indicies[i]);
            
            i = 0;
        }
        
        result.push_back(detail::generate_genotype(elements, element_indicies));
        ++element_indicies[i];
    }
    
    return result;
}

template <typename MappableType>
std::unordered_map<MappableType, unsigned> get_element_count_map(const Genotype<MappableType>& genotype)
{
    std::unordered_map<MappableType, unsigned> result {};
    result.reserve(genotype.zygosity());
    
    for (unsigned i {}; i < genotype.ploidy(); ++i) {
        ++result[genotype.at(i)];
    }
    
    return result;
}

template <typename MappableType>
std::ostream& operator<<(std::ostream& os, const Genotype<MappableType>& genotype)
{
    if (genotype.ploidy() == 0) {
        os << "empty genotype";
        return os;
    }
    auto element_counts = get_element_count_map(genotype);
    std::vector<std::pair<MappableType, unsigned>> p {element_counts.begin(), element_counts.end()};
    for (unsigned i {}; i < p.size() - 1; ++i) {
        os << p[i].first << "(" << p[i].second << "),";
    }
    os << p.back().first << "(" << p.back().second << ")";
    return os;
}

void print_alleles(const Genotype<Haplotype>& genotype);
void print_variant_alleles(const Genotype<Haplotype>& genotype);

#endif /* defined(__Octopus__genotype__) */
