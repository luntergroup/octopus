//
//  genotype.h
//  Octopus
//
//  Created by Daniel Cooke on 24/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__genotype__
#define __Octopus__genotype__

#include <vector>
#include <unordered_map>
#include <ostream>
#include <iterator>  // std::cbegin etc
#include <initializer_list>
#include <algorithm> // std::inplace_merge, std::all_of, std::binary_search, std::equal_range, std::unique_copy, std::equal
#include <boost/functional/hash.hpp> // boost::hash_range
#include <boost/math/special_functions/binomial.hpp>

#include "haplotype.h"
#include "equitable.h"

template <typename MappableType>
class Genotype : public Equitable<Genotype<MappableType>>
{
public:
    using Iterator = typename std::vector<MappableType>::const_iterator;
    
    Genotype() = default;
    Genotype(unsigned ploidy);
    Genotype(std::initializer_list<MappableType> the_elements);
    ~Genotype() = default;
    
    Genotype(const Genotype&)            = default;
    Genotype& operator=(const Genotype&) = default;
    Genotype(Genotype&&)                 = default;
    Genotype& operator=(Genotype&&)      = default;
    
    const MappableType& at(unsigned n) const;
    void emplace(const MappableType& an_element);
    void emplace(MappableType&& an_element);
    
    Iterator begin() const;
    Iterator end() const;
    Iterator cbegin() const;
    Iterator cend() const;
    
    unsigned ploidy() const noexcept;
    bool contains(const MappableType& an_element) const;
    unsigned num_occurences(const MappableType& an_element) const;
    bool is_homozygous() const;
    std::vector<MappableType> get_unique() const;
    
private:
    std::vector<MappableType> the_elements_;
};

template <typename MappableType>
Genotype<MappableType>::Genotype(unsigned ploidy)
:
the_elements_ {}
{
    the_elements_.reserve(ploidy);
}

template <typename MappableType>
Genotype<MappableType>::Genotype(std::initializer_list<MappableType> the_elements)
:
the_elements_ {the_elements}
{}

template <typename MappableType>
const MappableType& Genotype<MappableType>::at(unsigned n) const
{
    return the_elements_.at(n);
}

template <typename MappableType>
void Genotype<MappableType>::emplace(const MappableType& an_element)
{
    the_elements_.emplace_back(an_element);
    std::inplace_merge(std::begin(the_elements_), std::prev(std::end(the_elements_)), std::end(the_elements_));
}

template <typename MappableType>
void Genotype<MappableType>::emplace(MappableType&& an_element)
{
    the_elements_.emplace_back(std::move(an_element));
    std::inplace_merge(std::begin(the_elements_), std::prev(std::end(the_elements_)), std::end(the_elements_));
}

template <typename MappableType>
typename Genotype<MappableType>::Iterator Genotype<MappableType>::begin() const
{
    return the_elements_.begin();
}

template <typename MappableType>
typename Genotype<MappableType>::Iterator Genotype<MappableType>::end() const
{
    return the_elements_.end();
}

template <typename MappableType>
typename Genotype<MappableType>::Iterator Genotype<MappableType>::cbegin() const
{
    return the_elements_.cbegin();
}

template <typename MappableType>
typename Genotype<MappableType>::Iterator Genotype<MappableType>::cend() const
{
    return the_elements_.cend();
}

template <typename MappableType>
bool Genotype<MappableType>::is_homozygous() const
{
    const auto& first_element = the_elements_.front();
    return std::all_of(std::next(std::cbegin(the_elements_)), std::cend(the_elements_),
                       [&first_element] (const auto& an_element) {
                           return first_element == an_element;
                       });
}

template <typename MappableType>
bool Genotype<MappableType>::contains(const MappableType& an_element) const
{
    return std::binary_search(std::cbegin(the_elements_), std::cend(the_elements_), an_element);
}

template <typename MappableType>
unsigned Genotype<MappableType>::num_occurences(const MappableType& an_element) const
{
    auto equal_range = std::equal_range(std::cbegin(the_elements_), std::cend(the_elements_), an_element);
    return static_cast<unsigned>(std::distance(equal_range.first, equal_range.second));
}

template <typename MappableType>
unsigned Genotype<MappableType>::ploidy() const noexcept
{
    return static_cast<unsigned>(the_elements_.size());
}

template <typename MappableType>
std::vector<MappableType> Genotype<MappableType>::get_unique() const
{
    std::vector<MappableType> result {};
    result.reserve(ploidy());
    
    std::unique_copy(std::cbegin(the_elements_), std::cend(the_elements_), std::back_inserter(result));
    
    return result;
}

template <typename MappableType>
bool operator==(const Genotype<MappableType>& lhs, const Genotype<MappableType>& rhs)
{
    if (lhs.ploidy() != rhs.ploidy()) return false;
    return std::equal(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs));
}

namespace std {
    template <typename MappableType> struct hash<Genotype<MappableType>>
    {
        size_t operator()(const Genotype<MappableType>& genotype) const
        {
            return boost::hash_range(std::cbegin(genotype), std::cend(genotype));
        }
    };
}

template <typename MappableType1, typename MappableType2>
bool contains(const Genotype<MappableType1>& lhs, const Genotype<MappableType2>& rhs)
{
    return std::all_of(std::cbegin(lhs), std::cend(lhs), [&rhs] (const auto& element) {
        return rhs.contains(element);
    });
}

inline unsigned num_genotypes(unsigned num_elements, unsigned ploidy)
{
    return static_cast<unsigned>(boost::math::binomial_coefficient<double>(num_elements + ploidy - 1,
                                                                           num_elements - 1));
}

namespace detail
{
    template <typename MappableType>
    Genotype<MappableType> get_genotype_from_element_indicies(const std::vector<MappableType>& the_elements,
                                                              const std::vector<unsigned>& element_indicies)
    {
        Genotype<MappableType> result {};
        
        for (auto i : element_indicies) {
            result.emplace(the_elements.at(i));
        }
        
        return result;
    }
} // end namespace detail

// Assumes the input haplotypes are unique
template <typename MappableType>
std::vector<Genotype<MappableType>> get_all_genotypes(const std::vector<MappableType>& elements, unsigned ploidy)
{
    std::vector<Genotype<MappableType>> result {};
    
    unsigned num_elements {static_cast<unsigned>(elements.size())};
    
    result.reserve(num_genotypes(num_elements, ploidy));
    
    std::vector<unsigned> element_indicies(ploidy, 0);
    
    unsigned i {};
    while (true) {
        if (element_indicies[i] == num_elements) {
            do {
                ++i;
            } while (element_indicies[i] == num_elements - 1);
            
            if (i == ploidy) return result;
            
            ++element_indicies[i];
            
            for (unsigned j {0}; j <= i; ++j) {
                element_indicies[j] = element_indicies[i];
            }
            
            i = 0;
        }
        
        result.push_back(detail::get_genotype_from_element_indicies(elements, element_indicies));
        ++element_indicies[i];
    }
}

template <typename MappableType>
std::unordered_map<MappableType, unsigned> get_element_occurence_map(const Genotype<MappableType>& a_genotype)
{
    std::unordered_map<MappableType, unsigned> result {};
    
    for (unsigned i {}; i < a_genotype.ploidy(); ++i) {
        ++result[a_genotype.at(i)];
    }
    
    return result;
}

template <typename MappableType>
inline std::ostream& operator<<(std::ostream& os, const Genotype<MappableType>& a_genotype)
{
    auto element_occurences = get_element_occurence_map(a_genotype);
    std::vector<std::pair<Haplotype, unsigned>> p {element_occurences.begin(), element_occurences.end()};
    for (unsigned i {}; i < p.size() - 1; ++i) {
        os << p[i].first << "(" << p[i].second << "),";
    }
    os << p.back().first << "(" << p.back().second << ")";
    return os;
}

#endif /* defined(__Octopus__genotype__) */
