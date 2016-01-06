//
//  cancer_genotype.hpp
//  Octopus
//
//  Created by Daniel Cooke on 01/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef cancer_genotype_hpp
#define cancer_genotype_hpp

#include <initializer_list>
#include <utility>
#include <ostream>

#include <boost/functional/hash.hpp>

#include "equitable.hpp"
#include "genotype.hpp"

template <typename MappableType>
class CancerGenotype : public Equitable<CancerGenotype<MappableType>>
{
public:
    CancerGenotype() = default;
    explicit CancerGenotype(std::initializer_list<MappableType> normal_elements, MappableType cancer_element);
    explicit CancerGenotype(Genotype<MappableType> germline_genotype, MappableType cancer_element);
    ~CancerGenotype() = default;
    
    CancerGenotype(const CancerGenotype&)            = default;
    CancerGenotype& operator=(const CancerGenotype&) = default;
    CancerGenotype(CancerGenotype&&)                 = default;
    CancerGenotype& operator=(CancerGenotype&&)      = default;
    
    const MappableType& at(unsigned n) const;
    const MappableType& operator[](unsigned n) const;
    
    const Genotype<MappableType>& get_germline_genotype() const;
    const MappableType& get_cancer_element() const;
    
    unsigned ploidy() const noexcept;
    bool contains(const MappableType& element) const;
    unsigned count(const MappableType& element) const;
    bool is_homozygous() const;
    unsigned zygosity() const;
    std::vector<MappableType> get_unique() const;
    
private:
    Genotype<MappableType> germline_genotype_;
    MappableType cancer_element_;
};

template <typename MappableType>
CancerGenotype<MappableType>::CancerGenotype(std::initializer_list<MappableType> germline_elements, MappableType cancer_element)
:
germline_genotype_ {germline_elements},
cancer_element_ {std::move(cancer_element)}
{}

template <typename MappableType>
CancerGenotype<MappableType>::CancerGenotype(Genotype<MappableType> germline_genotype, MappableType cancer_element)
:
germline_genotype_ {std::move(germline_genotype)},
cancer_element_ {std::move(cancer_element)}
{}

template <typename MappableType>
const MappableType& CancerGenotype<MappableType>::at(unsigned n) const
{
    return (n < ploidy()) ? germline_genotype_.at(n) : cancer_element_;
}

template <typename MappableType>
const MappableType& CancerGenotype<MappableType>::operator[](unsigned n) const
{
    return (n < ploidy()) ? germline_genotype_[n] : cancer_element_;
}

template <typename MappableType>
const Genotype<MappableType>& CancerGenotype<MappableType>::get_germline_genotype() const
{
    return germline_genotype_;
}

template <typename MappableType>
const MappableType& CancerGenotype<MappableType>::get_cancer_element() const
{
    return cancer_element_;
}

template <typename MappableType>
unsigned CancerGenotype<MappableType>::ploidy() const noexcept
{
    return germline_genotype_.ploidy();
}

template <typename MappableType>
bool CancerGenotype<MappableType>::contains(const MappableType& element) const
{
    return germline_genotype_.contains(element) || cancer_element_ == element;
}

template <typename MappableType>
unsigned CancerGenotype<MappableType>::count(const MappableType& element) const
{
    return germline_genotype_.count(element) + ((cancer_element_ == element) ? 1 : 0);
}

template <typename MappableType>
bool CancerGenotype<MappableType>::is_homozygous() const
{
    return germline_genotype_.is_homozygous() && cancer_element_ == germline_genotype_[0];
}

template <typename MappableType>
unsigned CancerGenotype<MappableType>::zygosity() const
{
    return germline_genotype_.zygosity() + ((germline_genotype_.contains(cancer_element_)) ? 0 : 1);
}

template <typename MappableType>
std::vector<MappableType> CancerGenotype<MappableType>::get_unique() const
{
    auto result = germline_genotype_.get_unique();
    if (!germline_genotype_.contains(cancer_element_)) result.push_back(cancer_element_);
    return result;
}

// non-member methods

template <typename MappableType2, typename MappableType1>
CancerGenotype<MappableType2> splice(const CancerGenotype<MappableType1>& genotype, const GenomicRegion& region)
{
    return CancerGenotype<MappableType2> {
        splice<MappableType2>(genotype.get_germline_genotype(), region),
        splice<MappableType2>(genotype.get_cancer_element(), region)
    };
}

inline size_t num_cancer_genotypes(const unsigned num_elements, const unsigned ploidy)
{
    return (num_genotypes(num_elements, ploidy) - 1) * num_elements;
}

template <typename MappableType>
std::vector<CancerGenotype<MappableType>>
generate_all_cancer_genotypes(const std::vector<Genotype<MappableType>>& germline_genotypes,
                              const std::vector<MappableType>& elements)
{
    std::vector<CancerGenotype<MappableType>> result {};
    result.reserve((germline_genotypes.size() - 1) * elements.size());
    
    for (const auto& germline_genotype : germline_genotypes) {
        for (const auto& element : elements) {
            if (!is_homozygous(germline_genotype, element)) {
                result.emplace_back(germline_genotype, element);
            }
        }
    }
    
    return result;
}

template <typename MappableType>
std::vector<CancerGenotype<MappableType>>
generate_all_cancer_genotypes(const std::vector<MappableType>& elements, const unsigned ploidy)
{
    auto germline_genotypes = generate_all_genotypes(elements, ploidy);
    
    std::vector<CancerGenotype<MappableType>> result {};
    result.reserve(num_cancer_genotypes(static_cast<unsigned>(elements.size()), ploidy));
    
    for (const auto& germline_genotype : germline_genotypes) {
        for (const auto& element : elements) {
            if (!is_homozygous(germline_genotype, element)) {
                result.emplace_back(germline_genotype, element);
            }
        }
    }
    
    return result;
}

template <typename MappableType>
bool operator==(const CancerGenotype<MappableType>& lhs, const CancerGenotype<MappableType>& rhs)
{
    return lhs.get_cancer_element() == rhs.get_cancer_element() && lhs.get_germline_genotype() == rhs.get_germline_genotype();
}

namespace std
{
    template <typename MappableType> struct hash<CancerGenotype<MappableType>>
    {
        size_t operator()(const CancerGenotype<MappableType>& genotype) const
        {
            size_t seed {};
            boost::hash_combine(seed, hash<Genotype<MappableType>>()(genotype.get_germline_genotype()));
            boost::hash_combine(seed, hash<MappableType>()(genotype.get_cancer_element()));
            return seed;
        }
    };
} // namespace std

template <typename MappableType>
std::ostream& operator<<(std::ostream& os, const CancerGenotype<MappableType>& genotype)
{
    os << genotype.get_germline_genotype() << "," << genotype.get_cancer_element() << "(cancer)";
    return os;
}

inline void print_alleles(const CancerGenotype<Haplotype>& genotype)
{
    print_alleles(genotype.get_germline_genotype());
    std::cout << " + ";
    print_alleles(genotype.get_cancer_element());
}

inline void print_variant_alleles(const CancerGenotype<Haplotype>& genotype)
{
    print_variant_alleles(genotype.get_germline_genotype());
    std::cout << " + ";
    print_variant_alleles(genotype.get_cancer_element());
}

#endif /* cancer_genotype_h */
