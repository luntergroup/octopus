//
//  cancer_genotype.h
//  Octopus
//
//  Created by Daniel Cooke on 01/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef cancer_genotype_h
#define cancer_genotype_h

#include <initializer_list>
#include <utility>
#include <ostream>
#include <boost/functional/hash.hpp> // boost::hash_combine

#include "equitable.hpp"
#include "genotype.hpp"

template <typename MappableType>
class CancerGenotype : public Equitable<CancerGenotype<MappableType>>
{
public:
    CancerGenotype() = default;
    explicit CancerGenotype(std::initializer_list<MappableType> normal_elements, MappableType cancer_element);
    explicit CancerGenotype(Genotype<MappableType> normal_genotype, MappableType cancer_element);
    ~CancerGenotype() = default;
    
    CancerGenotype(const CancerGenotype&)            = default;
    CancerGenotype& operator=(const CancerGenotype&) = default;
    CancerGenotype(CancerGenotype&&)                 = default;
    CancerGenotype& operator=(CancerGenotype&&)      = default;
    
    const MappableType& at(unsigned n) const;
    const MappableType& operator[](unsigned n) const;
    
    const Genotype<MappableType>& get_normal_genotype() const;
    const MappableType& get_cancer_element() const;
    
    unsigned ploidy() const noexcept;
    bool contains(const MappableType& element) const;
    unsigned num_occurences(const MappableType& element) const;
    bool is_homozygous() const;
    unsigned zygosity() const;
    std::vector<MappableType> get_unique() const;
    
private:
    Genotype<MappableType> normal_genotype_;
    MappableType cancer_element_;
};

template <typename MappableType>
CancerGenotype<MappableType>::CancerGenotype(std::initializer_list<MappableType> elements, MappableType cancer_element)
:
normal_genotype_ {elements},
cancer_element_ {std::move(cancer_element)}
{}

template <typename MappableType>
CancerGenotype<MappableType>::CancerGenotype(Genotype<MappableType> normal_genotype, MappableType cancer_element)
:
normal_genotype_ {std::move(normal_genotype)},
cancer_element_ {std::move(cancer_element)}
{}

template <typename MappableType>
const MappableType& CancerGenotype<MappableType>::at(unsigned n) const
{
    return (n < ploidy()) ? normal_genotype_.at(n) : cancer_element_;
}

template <typename MappableType>
const MappableType& CancerGenotype<MappableType>::operator[](unsigned n) const
{
    return (n < ploidy()) ? normal_genotype_[n] : cancer_element_;
}

template <typename MappableType>
const Genotype<MappableType>& CancerGenotype<MappableType>::get_normal_genotype() const
{
    return normal_genotype_;
}

template <typename MappableType>
const MappableType& CancerGenotype<MappableType>::get_cancer_element() const
{
    return cancer_element_;
}

template <typename MappableType>
unsigned CancerGenotype<MappableType>::ploidy() const noexcept
{
    return normal_genotype_.ploidy();
}

template <typename MappableType>
bool CancerGenotype<MappableType>::contains(const MappableType& element) const
{
    return normal_genotype_.contains(element) || cancer_element_ == element;
}

template <typename MappableType>
unsigned CancerGenotype<MappableType>::num_occurences(const MappableType& element) const
{
    return normal_genotype_.num_occurences(element) + (cancer_element_ == element) ? 1 : 0;
}

template <typename MappableType>
bool CancerGenotype<MappableType>::is_homozygous() const
{
    return normal_genotype_.is_homozygous() && cancer_element_ == normal_genotype_[0];
}

template <typename MappableType>
unsigned CancerGenotype<MappableType>::zygosity() const
{
    return normal_genotype_.zygosity() + (normal_genotype_.contains(cancer_element_)) ? 0 : 1;
}

template <typename MappableType>
std::vector<MappableType> CancerGenotype<MappableType>::get_unique() const
{
    auto result = normal_genotype_.get_unique();
    if (!normal_genotype_.contains(cancer_element_)) result.push_back(cancer_element_);
    return result;
}

template <typename MappableType>
bool operator==(const CancerGenotype<MappableType>& lhs, const CancerGenotype<MappableType>& rhs)
{
    return lhs.get_cancer_element() == rhs.get_cancer_element() && lhs.get_normal_genotype() == rhs.get_normal_genotype();
}

namespace std
{
    template <typename MappableType> struct hash<CancerGenotype<MappableType>>
    {
        size_t operator()(const CancerGenotype<MappableType>& genotype) const
        {
            size_t seed {};
            boost::hash_combine(seed, hash<Genotype<MappableType>>()(genotype.get_normal_genotype()));
            boost::hash_combine(seed, hash<MappableType>()(genotype.get_cancer_element()));
            return seed;
        }
    };
} // namespace std

template <typename MappableType>
std::vector<CancerGenotype<MappableType>>
generate_all_cancer_genotypes(const std::vector<MappableType>& elements, unsigned ploidy)
{
    auto normal_genotypes = generate_all_genotypes(elements, ploidy);
    
    std::vector<CancerGenotype<MappableType>> result {};
    result.reserve(normal_genotypes.size() * (elements.size() - 1));
    
    for (auto normal_genotype : normal_genotypes) {
        for (const auto& element : elements) {
            if (!normal_genotype.contains(element)) {
                result.emplace_back(normal_genotype, element);
            }
        }
    }
    
    return result;
}

template <typename MappableType>
std::ostream& operator<<(std::ostream& os, const CancerGenotype<MappableType>& genotype)
{
    os << genotype.get_normal_genotype() << "," << genotype.get_cancer_element() << "(cancer)";
    return os;
}

inline void print_alleles(const CancerGenotype<Haplotype>& genotype)
{
    print_alleles(genotype.get_normal_genotype());
    std::cout << std::endl;
    print_alleles(genotype.get_cancer_element());
}

#endif /* cancer_genotype_h */
