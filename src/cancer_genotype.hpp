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
#include <memory>

#include <boost/functional/hash.hpp>

#include "equitable.hpp"
#include "genotype.hpp"

namespace Octopus
{
template <typename MappableType>
class CancerGenotype : public Equitable<CancerGenotype<MappableType>>
{
public:
    CancerGenotype() = default;
    
    explicit CancerGenotype(std::initializer_list<MappableType> normal_elements,
                            const MappableType& cancer_element);
    explicit CancerGenotype(std::initializer_list<MappableType> normal_elements,
                            MappableType&& cancer_element);
    
    template <typename G>
    explicit CancerGenotype(G&& germline_genotype, const std::shared_ptr<MappableType>& cancer_element);
    
    template <typename G, typename C>
    explicit CancerGenotype(G&& germline_genotype, C&& cancer_element);
    
    ~CancerGenotype() = default;
    
    CancerGenotype(const CancerGenotype&)            = default;
    CancerGenotype& operator=(const CancerGenotype&) = default;
    CancerGenotype(CancerGenotype&&)                 = default;
    CancerGenotype& operator=(CancerGenotype&&)      = default;
    
    const MappableType& operator[](unsigned n) const;
    
    const Genotype<MappableType>& get_germline_genotype() const;
    const MappableType& get_cancer_element() const;
    
    unsigned ploidy() const noexcept;
    bool contains(const MappableType& element) const;
    unsigned count(const MappableType& element) const;
    bool is_homozygous() const;
    unsigned zygosity() const;
    std::vector<MappableType> copy_unique() const;
    
private:
    Genotype<MappableType> germline_genotype_;
    std::shared_ptr<MappableType> cancer_element_;
};

template <typename MappableType>
CancerGenotype<MappableType>::CancerGenotype(std::initializer_list<MappableType> germline_elements,
                                             const MappableType& cancer_element)
:
germline_genotype_ {germline_elements},
cancer_element_ {std::make_shared<MappableType>(cancer_element)}
{}

template <typename MappableType>
CancerGenotype<MappableType>::CancerGenotype(std::initializer_list<MappableType> germline_elements,
                                             MappableType&& cancer_element)
:
germline_genotype_ {germline_elements},
cancer_element_ {std::make_shared<MappableType>(std::move(cancer_element))}
{}

template <typename MappableType>
template <typename G>
CancerGenotype<MappableType>::CancerGenotype(G&& germline_genotype,
                                             const std::shared_ptr<MappableType>& cancer_element)
:
germline_genotype_ {std::forward<G>(germline_genotype)},
cancer_element_ {cancer_element}
{}

template <typename MappableType>
template <typename G, typename C>
CancerGenotype<MappableType>::CancerGenotype(G&& germline_genotype, C&& cancer_element)
:
germline_genotype_ {std::forward<G>(germline_genotype)},
cancer_element_ {std::make_shared<MappableType>(std::forward<C>(cancer_element))}
{}

template <typename MappableType>
const MappableType& CancerGenotype<MappableType>::operator[](unsigned n) const
{
    return (n < ploidy()) ? germline_genotype_[n] : *cancer_element_;
}

template <typename MappableType>
const Genotype<MappableType>& CancerGenotype<MappableType>::get_germline_genotype() const
{
    return germline_genotype_;
}

template <typename MappableType>
const MappableType& CancerGenotype<MappableType>::get_cancer_element() const
{
    return *cancer_element_;
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
    return germline_genotype_.count(element) + ((*cancer_element_ == element) ? 1 : 0);
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
std::vector<MappableType> CancerGenotype<MappableType>::copy_unique() const
{
    auto result = germline_genotype_.get_unique();
    if (!germline_genotype_.contains(cancer_element_)) result.push_back(cancer_element_);
    return result;
}

// non-member methods

template <typename MappableType2, typename MappableType1>
CancerGenotype<MappableType2> splice(const CancerGenotype<MappableType1>& genotype,
                                     const GenomicRegion& region)
{
    return CancerGenotype<MappableType2> {
        splice<MappableType2>(genotype.get_germline_genotype(), region),
        splice<MappableType2>(genotype.get_cancer_element(), region)
    };
}

std::pair<std::vector<CancerGenotype<Haplotype>>, std::vector<Genotype<Haplotype>>>
generate_all_cancer_genotypes(const std::vector<Haplotype>& haplotypes, const unsigned ploidy);

template <typename MappableType>
bool operator==(const CancerGenotype<MappableType>& lhs, const CancerGenotype<MappableType>& rhs)
{
    return lhs.get_cancer_element() == rhs.get_cancer_element()
                && lhs.get_germline_genotype() == rhs.get_germline_genotype();
}
} // namespace Octopus

namespace std
{
    template <typename MappableType> struct hash<Octopus::CancerGenotype<MappableType>>
    {
        size_t operator()(const Octopus::CancerGenotype<MappableType>& genotype) const
        {
            using boost::hash_combine;
            size_t result {0};
            hash_combine(result, hash<Genotype<MappableType>>()(genotype.get_germline_genotype()));
            hash_combine(result, hash<MappableType>()(genotype.get_cancer_element()));
            return result;
        }
    };
    
    template <typename MappableType>
    struct hash<reference_wrapper<const Octopus::CancerGenotype<MappableType>>>
    {
        size_t operator()(const reference_wrapper<const Octopus::CancerGenotype<MappableType>> genotype) const
        {
            return hash<Octopus::CancerGenotype<MappableType>>()(genotype);
        }
    };
} // namespace std

namespace Octopus
{
template <typename MappableType>
std::ostream& operator<<(std::ostream& os, const CancerGenotype<MappableType>& genotype)
{
    os << genotype.get_germline_genotype() << "," << genotype.get_cancer_element() << "(cancer)";
    return os;
}

namespace debug
{
    template <typename S>
    void print_alleles(S&& stream, const CancerGenotype<Haplotype>& genotype)
    {
        ::debug::print_alleles(stream, genotype.get_germline_genotype());
        stream << " + ";
        ::debug::print_alleles(stream, genotype.get_cancer_element());
    }
    
    void print_alleles(const CancerGenotype<Haplotype>& genotype);
    
    template <typename S>
    void print_variant_alleles(S&& stream, const CancerGenotype<Haplotype>& genotype)
    {
        ::debug::print_variant_alleles(stream, genotype.get_germline_genotype());
        stream << " + ";
        ::debug::print_variant_alleles(stream, genotype.get_cancer_element());
    }
    
    void print_variant_alleles(const CancerGenotype<Haplotype>& genotype);
} // namespace debug
} // namespace Octopus

#endif /* cancer_genotype_h */
