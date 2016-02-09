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

#include <boost/functional/hash.hpp>

#include "allele.hpp"
#include "haplotype.hpp"
#include "equitable.hpp"
#include "mappable.hpp"

template <typename T>
using EnableIfGenotypable = std::enable_if_t<std::is_same<T, Haplotype>::value || std::is_same<T, Allele>::value>;

template <typename MappableType, typename = EnableIfGenotypable<MappableType>>
class Genotype;

template <typename MappableType>
class Genotype<MappableType> : public Equitable<Genotype<MappableType>>, public Mappable<Genotype<MappableType>>
{
public:
    Genotype() = default;
    
    explicit Genotype(unsigned ploidy);
    explicit Genotype(unsigned ploidy, const MappableType& init);
    explicit Genotype(unsigned ploidy, const std::shared_ptr<MappableType>& init);
    explicit Genotype(std::initializer_list<MappableType> elements);
    explicit Genotype(std::initializer_list<std::shared_ptr<MappableType>> elements);
    
    ~Genotype() = default;
    
    Genotype(const Genotype&)            = default;
    Genotype& operator=(const Genotype&) = default;
    Genotype(Genotype&&)                 = default;
    Genotype& operator=(Genotype&&)      = default;
    
    template <typename T> void emplace(T&& element);
    void emplace(const std::shared_ptr<MappableType>& element);
    
    const MappableType& at(unsigned n) const;
    const MappableType& operator[](unsigned n) const;
    
    const GenomicRegion& get_region() const noexcept;
    
    unsigned ploidy() const noexcept;
    
    bool contains(const MappableType& element) const;
    unsigned count(const MappableType& element) const;
    bool is_homozygous() const;
    unsigned zygosity() const;
    
    std::vector<MappableType> copy_unique() const;
    
private:
    using ElementPtr    = std::shared_ptr<MappableType>;
    using BaseContainer = std::vector<ElementPtr>;
    using BaseIterator  = typename BaseContainer::const_iterator;
    
    BaseContainer elements_;
    
    struct ElementLess
    {
        bool operator()(const ElementPtr& lhs, const ElementPtr& rhs) const
        {
            return *lhs < *rhs;
        }
        bool operator()(const MappableType& lhs, const ElementPtr& rhs) const
        {
            return lhs < *rhs;
        }
        bool operator()(const ElementPtr& lhs, const MappableType& rhs) const
        {
            return *lhs < rhs;
        }
    };
    
    struct ElementEqual
    {
        bool operator()(const ElementPtr& lhs, const ElementPtr& rhs) const
        {
            return *lhs == *rhs;
        }
        bool operator()(const MappableType& lhs, const ElementPtr& rhs) const
        {
            return lhs == *rhs;
        }
        bool operator()(const ElementPtr& lhs, const MappableType& rhs) const
        {
            return *lhs == rhs;
        }
    };
    
public:
    class Iterator : public BaseIterator
    {
    public:
        using value_type = MappableType;
        using reference  = const MappableType&;
        using pointer    = const MappableType*;
        
        Iterator(BaseIterator it) : BaseIterator {it} {}
        
        reference operator*() const
        {
            return *BaseIterator::operator*();
        }
    };
    
    Iterator begin() const noexcept;
    Iterator end() const noexcept;
    Iterator cbegin() const noexcept;
    Iterator cend() const noexcept;
};

template <>
class Genotype<Allele> : public Equitable<Genotype<Allele>>, public Mappable<Genotype<Allele>>
{
public:
    Genotype() = default;
    explicit Genotype(unsigned ploidy);
    explicit Genotype(unsigned ploidy, const Allele& init);
    explicit Genotype(std::initializer_list<Allele> alleles);
    ~Genotype() = default;
    
    Genotype(const Genotype&)            = default;
    Genotype& operator=(const Genotype&) = default;
    Genotype(Genotype&&)                 = default;
    Genotype& operator=(Genotype&&)      = default;
    
    template <typename T> void emplace(T&& allele);
    
    const GenomicRegion& get_region() const noexcept;
    
    const Allele& at(unsigned n) const;
    const Allele& operator[](unsigned n) const;
    
    unsigned ploidy() const noexcept;
    
    bool contains(const Allele& allele) const;
    unsigned count(const Allele& allele) const;
    bool is_homozygous() const;
    unsigned zygosity() const;
    
    std::vector<Allele> copy_unique() const;
    
private:
    std::vector<Allele> alleles_;
    
public:
    using Iterator = decltype(alleles_)::const_iterator;
    
    Iterator begin() const noexcept;
    Iterator end() const noexcept ;
    Iterator cbegin() const noexcept ;
    Iterator cend() const noexcept ;
};

// Genotype<MappableType>

template <typename MappableType>
Genotype<MappableType>::Genotype(const unsigned ploidy)
:
elements_ {}
{
    elements_.reserve(ploidy);
}

template <typename MappableType>
Genotype<MappableType>::Genotype(const unsigned ploidy, const MappableType& init)
{
    if (ploidy > 0) {
        elements_.resize(ploidy);
        elements_.front() = std::make_shared<MappableType>(init);
        std::fill_n(std::next(std::begin(elements_)), ploidy - 1, elements_.front());
    }
}

template <typename MappableType>
Genotype<MappableType>::Genotype(const unsigned ploidy, const std::shared_ptr<MappableType>& init)
:
elements_ {ploidy, init}
{}

template <typename MappableType>
Genotype<MappableType>::Genotype(std::initializer_list<MappableType> elements)
{
    std::transform(std::cbegin(elements), std::cend(elements), std::back_inserter(elements_),
                   [] (const auto& element) { return std::make_shared<MappableType>(element); });
    std::sort(std::begin(elements_), std::end(elements_), ElementLess {});
}

template <typename MappableType>
Genotype<MappableType>::Genotype(std::initializer_list<std::shared_ptr<MappableType>> elements)
:
elements_ {elements}
{
    std::sort(std::begin(elements_), std::end(elements_), ElementLess {});
}

template <typename MappableType>
template <typename T>
void Genotype<MappableType>::emplace(T&& element)
{
    elements_.emplace_back(std::make_shared<MappableType>(std::forward<T>(element)));
    std::inplace_merge(std::begin(elements_), std::prev(std::end(elements_)), std::end(elements_),
                       ElementLess {});
}

template <typename MappableType>
void Genotype<MappableType>::emplace(const std::shared_ptr<MappableType>& element)
{
    elements_.emplace_back(element);
    std::inplace_merge(std::begin(elements_), std::prev(std::end(elements_)), std::end(elements_),
                       ElementLess {});
}

template <typename MappableType>
const MappableType& Genotype<MappableType>::at(const unsigned n) const
{
    return *elements_.at(n);
}

template <typename MappableType>
const MappableType& Genotype<MappableType>::operator[](const unsigned n) const
{
    return *elements_[n];
}

template <typename MappableType>
const GenomicRegion& Genotype<MappableType>::get_region() const noexcept
{
    return elements_.front()->get_region();
}

template <typename MappableType>
unsigned Genotype<MappableType>::ploidy() const noexcept
{
    return static_cast<unsigned>(elements_.size());
}

template <typename MappableType>
bool Genotype<MappableType>::is_homozygous() const
{
    return *elements_.front() == *elements_.back();
}

template <typename MappableType>
unsigned Genotype<MappableType>::zygosity() const
{
    unsigned result {0};
    
    for (auto it = std::cbegin(elements_), last = std::cend(elements_); it != last; ++result) {
        // naive algorithm faster in practice than binary searching
        it = std::find_if_not(std::next(it), last, [it] (const auto& x) { return *x == **it; });
    }
    
    return result;
}

template <typename MappableType>
bool Genotype<MappableType>::contains(const MappableType& element) const
{
    return std::binary_search(std::cbegin(elements_), std::cend(elements_), element,
                              ElementLess {});
}

template <typename MappableType>
unsigned Genotype<MappableType>::count(const MappableType& element) const
{
    const auto equal_range = std::equal_range(std::cbegin(elements_), std::cend(elements_), element,
                                              ElementLess {});
    return static_cast<unsigned>(std::distance(equal_range.first, equal_range.second));
}

template <typename MappableType>
std::vector<MappableType> Genotype<MappableType>::copy_unique() const
{
    std::vector<std::reference_wrapper<const ElementPtr>> ptr_copy {};
    ptr_copy.reserve(ploidy());
    
    std::unique_copy(std::cbegin(elements_), std::cend(elements_), std::back_inserter(ptr_copy),
                     ElementEqual {});
    
    std::vector<MappableType> result {};
    result.reserve(ptr_copy.size());
    
    std::transform(std::cbegin(ptr_copy), std::cend(ptr_copy), std::back_inserter(result),
                   [] (const auto& ptr) { return *ptr.get(); });
    
    return result;
}

template <typename MappableType>
typename Genotype<MappableType>::Iterator Genotype<MappableType>::begin() const noexcept
{
    return Iterator {std::cbegin(elements_)};
}

template <typename MappableType>
typename Genotype<MappableType>::Iterator Genotype<MappableType>::end() const noexcept
{
    return Iterator {std::cend(elements_)};
}

template <typename MappableType>
typename Genotype<MappableType>::Iterator Genotype<MappableType>::cbegin() const noexcept
{
    return Iterator {std::cbegin(elements_)};
}

template <typename MappableType>
typename Genotype<MappableType>::Iterator Genotype<MappableType>::cend() const noexcept
{
    return Iterator {std::cend(elements_)};
}

// Genotype<Allele>

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
    return splice<MappableType2>(lhs, rhs.get_region()) == rhs;
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

struct GenotypeSortCompare
{
    template <typename T>
    bool operator()(const Genotype<T>& lhs, const Genotype<T>& rhs) const
    {
        return std::lexicographical_compare(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs), std::cend(rhs));
    }
};

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
            GenotypeTp {elements.front(), elements.front()},
            GenotypeTp {elements.front(), elements.back()},
            GenotypeTp {elements.back(), elements.back()}
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
    auto generate_genotype(const Container& elements,
                           const std::vector<unsigned>& element_indicies)
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
        
        if (ploidy == 0 || elements.empty()) return ResultType {};
        
        const auto num_elements = static_cast<unsigned>(elements.size());
        
        if (num_elements == 1) return ResultType {GenotypeTp {ploidy, elements.front()}};
        
        if (ploidy == 2) {
            if (num_elements == 2) return detail::generate_all_diploid_biallelic_genotypes(elements);
            if (num_elements == 3) return detail::generate_all_diploid_triallelic_genotypes(elements);
            if (num_elements == 4) return detail::generate_all_diploid_tetraallelic_genotypes(elements);
        }
        
        if (ploidy == 1) return detail::generate_all_haploid_genotypes(elements);
        
        ResultType result {};
        result.reserve(num_genotypes(num_elements, ploidy));
        
        std::vector<unsigned> element_indicies(ploidy, 0);
        
        unsigned i {0};
        
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
    struct ShouldShareMemory : public std::is_same<MappableType, Haplotype> {};
    
    template <typename MappableType>
    auto generate_all_genotypes(const std::vector<MappableType>& elements, const unsigned ploidy, std::true_type)
    {
        std::vector<std::shared_ptr<MappableType>> temp_pointers(elements.size());
        std::transform(std::cbegin(elements), std::cend(elements), std::begin(temp_pointers),
                       [] (const auto& element) { return std::make_shared<MappableType>(element); });
        return do_generate_all_genotypes(temp_pointers, ploidy);
    }
    
    template <typename MappableType>
    auto generate_all_genotypes(const std::vector<MappableType>& elements, const unsigned ploidy, std::false_type)
    {
        return do_generate_all_genotypes(elements, ploidy);
    }
} // namespace detail

template <typename MappableType>
auto generate_all_genotypes(const std::vector<MappableType>& elements, const unsigned ploidy)
{
    return detail::generate_all_genotypes(elements, ploidy, detail::ShouldShareMemory<MappableType> {});
}

namespace detail
{
    inline size_t estimate_num_elements(const size_t num_genotypes)
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
auto extract_all_elements(const std::vector<Genotype<MappableType>>& genotypes, const GenomicRegion)
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
auto get_element_count_map(const Genotype<MappableType>& genotype)
{
    std::unordered_map<MappableType, unsigned> result {};
    result.reserve(genotype.zygosity());
    
    for (unsigned i {0}; i < genotype.ploidy(); ++i) {
        ++result[genotype.at(i)];
    }
    
    return result;
}

template <typename MappableType>
auto get_element_ref_count_map(const Genotype<MappableType>& genotype)
{
    std::unordered_map<std::reference_wrapper<const MappableType>, unsigned> result {};
    result.reserve(genotype.zygosity());
    
    for (unsigned i {0}; i < genotype.ploidy(); ++i) {
        ++result[genotype.at(i)];
    }
    
    return result;
}

template <typename MappableType2, typename MappableType1>
auto splice_all(const std::vector<Genotype<MappableType1>>& genotypes, const GenomicRegion& region)
{
    std::vector<Genotype<MappableType2>> result {};
    result.reserve(genotypes.size());
    
    for (const auto& genotype : genotypes) {
        result.emplace_back(splice<MappableType2>(genotype, region));
    }
    
    std::sort(std::begin(result), std::end(result), GenotypeSortCompare {});
    
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
