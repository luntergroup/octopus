// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef cancer_genotype_hpp
#define cancer_genotype_hpp

#include <initializer_list>
#include <utility>
#include <memory>
#include <functional>

#include <boost/functional/hash.hpp>

#include <concepts/equitable.hpp>

#include "genotype.hpp"

namespace octopus {

template <typename MappableType>
class CancerGenotype
    : public Equitable<CancerGenotype<MappableType>>, public Mappable<CancerGenotype<MappableType>>
{
public:
    using MappingDomain = typename Genotype<MappableType>::MappingDomain;
    
    CancerGenotype() = default;
    
    CancerGenotype(std::initializer_list<MappableType> normal_elements,
                   const MappableType& somatic_element);
    CancerGenotype(std::initializer_list<MappableType> normal_elements,
                   MappableType&& somatic_element);
    
    template <typename G>
    CancerGenotype(G&& germline_genotype, const std::shared_ptr<MappableType>& somatic_element);
    
    template <typename G, typename C>
    CancerGenotype(G&& germline_genotype, C&& somatic_element);
    
    CancerGenotype(const CancerGenotype&)            = default;
    CancerGenotype& operator=(const CancerGenotype&) = default;
    CancerGenotype(CancerGenotype&&)                 = default;
    CancerGenotype& operator=(CancerGenotype&&)      = default;
    
    ~CancerGenotype() = default;
    
    const GenomicRegion& mapped_region() const noexcept;
    
    const MappableType& operator[](unsigned n) const;
    
    const Genotype<MappableType>& germline_genotype() const;
    
    const MappableType& somatic_element() const;
    
    unsigned ploidy() const noexcept;
    
    bool contains(const MappableType& element) const;
    
    unsigned count(const MappableType& element) const;
    
    bool is_homozygous() const;
    
    unsigned zygosity() const;
    
    std::vector<MappableType> copy_unique() const;
    
private:
    Genotype<MappableType> germline_genotype_;
    std::shared_ptr<MappableType> somatic_element_;
};

template <typename MappableType>
CancerGenotype<MappableType>::CancerGenotype(std::initializer_list<MappableType> germline_elements,
                                             const MappableType& somatic_element)
: germline_genotype_ {germline_elements}
, somatic_element_ {std::make_shared<MappableType>(somatic_element)}
{}

template <typename MappableType>
CancerGenotype<MappableType>::CancerGenotype(std::initializer_list<MappableType> germline_elements,
                                             MappableType&& somatic_element)
: germline_genotype_ {germline_elements}
, somatic_element_ {std::make_shared<MappableType>(std::move(somatic_element))}
{}

template <typename MappableType>
template <typename G>
CancerGenotype<MappableType>::CancerGenotype(G&& germline_genotype,
                                             const std::shared_ptr<MappableType>& somatic_element)
: germline_genotype_ {std::forward<G>(germline_genotype)}
, somatic_element_ {somatic_element}
{}

template <typename MappableType>
template <typename G, typename C>
CancerGenotype<MappableType>::CancerGenotype(G&& germline_genotype, C&& somatic_element)
: germline_genotype_ {std::forward<G>(germline_genotype)}
, somatic_element_ {std::make_shared<MappableType>(std::forward<C>(somatic_element))}
{}

template <typename MappableType>
const GenomicRegion& CancerGenotype<MappableType>::mapped_region() const noexcept
{
    using octopus::mapped_region;
    return mapped_region(*somatic_element_);
}

template <typename MappableType>
const MappableType& CancerGenotype<MappableType>::operator[](unsigned n) const
{
    return (n < ploidy()) ? germline_genotype_[n] : *somatic_element_;
}

template <typename MappableType>
const Genotype<MappableType>& CancerGenotype<MappableType>::germline_genotype() const
{
    return germline_genotype_;
}

template <typename MappableType>
const MappableType& CancerGenotype<MappableType>::somatic_element() const
{
    return *somatic_element_;
}

template <typename MappableType>
unsigned CancerGenotype<MappableType>::ploidy() const noexcept
{
    return germline_genotype_.ploidy();
}

template <typename MappableType>
bool CancerGenotype<MappableType>::contains(const MappableType& element) const
{
    return germline_genotype_.contains(element) || somatic_element_ == element;
}

template <typename MappableType>
unsigned CancerGenotype<MappableType>::count(const MappableType& element) const
{
    return germline_genotype_.count(element) + ((*somatic_element_ == element) ? 1 : 0);
}

template <typename MappableType>
bool CancerGenotype<MappableType>::is_homozygous() const
{
    return germline_genotype_.is_homozygous() && somatic_element_ == germline_genotype_[0];
}

template <typename MappableType>
unsigned CancerGenotype<MappableType>::zygosity() const
{
    return germline_genotype_.zygosity() + ((germline_genotype_.contains(somatic_element_)) ? 0 : 1);
}

template <typename MappableType>
std::vector<MappableType> CancerGenotype<MappableType>::copy_unique() const
{
    auto result = germline_genotype_.get_unique();
    if (!germline_genotype_.contains(somatic_element_)) result.push_back(somatic_element_);
    return result;
}

// non-member methods

bool contains(const CancerGenotype<Haplotype>& genotype, const Allele& allele);
bool includes(const CancerGenotype<Haplotype>& genotype, const Allele& allele);

template <typename MappableType2, typename MappableType1>
CancerGenotype<MappableType2> splice(const CancerGenotype<MappableType1>& genotype,
                                     const GenomicRegion& region)
{
    return CancerGenotype<MappableType2> {
        splice<MappableType2>(genotype.germline_genotype(), region),
        splice<MappableType2>(genotype.somatic_element(), region)
    };
}

template <typename MappableType1, typename MappableType2>
bool contains(const CancerGenotype<MappableType1>& lhs, const CancerGenotype<MappableType2>& rhs)
{
    return splice<MappableType2>(lhs, rhs.mapped_region()) == rhs;
}

std::pair<std::vector<CancerGenotype<Haplotype>>, std::vector<Genotype<Haplotype>>>
generate_all_cancer_genotypes(const std::vector<Haplotype>& haplotypes, const unsigned ploidy);

template <typename MappableType>
Genotype<MappableType> convert(const CancerGenotype<MappableType>& genotype)
{
    Genotype<MappableType> result {genotype.ploidy() + 1};
    
    for (const auto& e : genotype.germline_genotype()) {
        result.emplace(e);
    }
    
    result.emplace(genotype.somatic_element());
    
    return result;
}
    
template <typename MappableType>
bool operator==(const CancerGenotype<MappableType>& lhs, const CancerGenotype<MappableType>& rhs)
{
    return lhs.somatic_element() == rhs.somatic_element()
                && lhs.germline_genotype() == rhs.germline_genotype();
}

template <typename MappableType>
std::ostream& operator<<(std::ostream& os, const CancerGenotype<MappableType>& genotype)
{
    os << genotype.germline_genotype() << "," << genotype.somatic_element() << "(cancer)";
    return os;
}

struct CancerGenotypeHash
{
    template <typename T>
    std::size_t operator()(const CancerGenotype<T>& genotype) const
    {
        using boost::hash_combine;
        size_t result {};
        hash_combine(result, std::hash<Genotype<T>>()(genotype.germline_genotype()));
        hash_combine(result, std::hash<T>()(genotype.somatic_element()));
        return result;
    }
};

namespace debug {
    template <typename S>
    void print_alleles(S&& stream, const CancerGenotype<Haplotype>& genotype)
    {
        print_alleles(stream, genotype.germline_genotype());
        stream << " + ";
        print_alleles(stream, genotype.somatic_element());
    }
    
    void print_alleles(const CancerGenotype<Haplotype>& genotype);
    
    template <typename S>
    void print_variant_alleles(S&& stream, const CancerGenotype<Haplotype>& genotype)
    {
        print_variant_alleles(stream, genotype.germline_genotype());
        stream << " + ";
        print_variant_alleles(stream, genotype.somatic_element());
    }
    
    void print_variant_alleles(const CancerGenotype<Haplotype>& genotype);
} // namespace debug
} // namespace octopus

namespace std {
    template <typename MappableType> struct hash<octopus::CancerGenotype<MappableType>>
    {
        size_t operator()(const octopus::CancerGenotype<MappableType>& genotype) const
        {
            return octopus::CancerGenotypeHash()(genotype);
        }
    };
    
    template <typename MappableType>
    struct hash<reference_wrapper<const octopus::CancerGenotype<MappableType>>>
    {
        size_t operator()(const reference_wrapper<const octopus::CancerGenotype<MappableType>> genotype) const
        {
            return hash<octopus::CancerGenotype<MappableType>>()(genotype);
        }
    };
} // namespace std

#endif
