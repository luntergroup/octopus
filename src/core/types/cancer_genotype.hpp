// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef cancer_genotype_hpp
#define cancer_genotype_hpp

#include <initializer_list>
#include <utility>
#include <memory>
#include <functional>
#include <iterator>
#include <algorithm>

#include <boost/functional/hash.hpp>

#include "concepts/equitable.hpp"
#include "utils/append.hpp"
#include "genotype.hpp"

namespace octopus {

template <typename MappableType>
class CancerGenotype
    : public Equitable<CancerGenotype<MappableType>>, public Mappable<CancerGenotype<MappableType>>
{
public:
    using MappingDomain = typename Genotype<MappableType>::MappingDomain;
    
    CancerGenotype() = default;
    
    CancerGenotype(std::initializer_list<MappableType> germline,
                   std::initializer_list<MappableType> somatic);
    CancerGenotype(std::initializer_list<MappableType> germline,
                   const MappableType& somatic);
    CancerGenotype(std::initializer_list<MappableType> germline,
                   MappableType&& somatic);
    template <typename G, typename S>
    CancerGenotype(G&& germline, S&& somatic);
    template <typename G>
    CancerGenotype(G&& germline, const std::shared_ptr<MappableType>& somatic);
    
    CancerGenotype(const CancerGenotype&)            = default;
    CancerGenotype& operator=(const CancerGenotype&) = default;
    CancerGenotype(CancerGenotype&&)                 = default;
    CancerGenotype& operator=(CancerGenotype&&)      = default;
    
    ~CancerGenotype() = default;
    
    const GenomicRegion& mapped_region() const noexcept;
    
    const MappableType& operator[](unsigned n) const;
    
    const Genotype<MappableType>& germline() const;
    const Genotype<MappableType>& somatic() const;
    
    unsigned germline_ploidy() const noexcept;
    unsigned somatic_ploidy() const noexcept;
    unsigned ploidy() const noexcept;
    bool contains(const MappableType& element) const;
    unsigned count(const MappableType& element) const;
    bool is_homozygous() const;
    unsigned zygosity() const;
    
    std::vector<MappableType> copy_unique() const;
    
private:
    Genotype<MappableType> germline_, somatic_;
};

template <typename MappableType>
CancerGenotype<MappableType>::CancerGenotype(std::initializer_list<MappableType> germline, std::initializer_list<MappableType> somatic)
: germline_ {germline}
, somatic_ {somatic}
{}

template <typename MappableType>
CancerGenotype<MappableType>::CancerGenotype(std::initializer_list<MappableType> germline, const MappableType& somatic)
: germline_ {germline}
, somatic_ {std::make_shared<MappableType>(somatic)}
{}

template <typename MappableType>
CancerGenotype<MappableType>::CancerGenotype(std::initializer_list<MappableType> germline, MappableType&& somatic)
: germline_ {germline}
, somatic_ {std::make_shared<MappableType>(std::move(somatic))}
{}

template <typename MappableType>
template <typename G>
CancerGenotype<MappableType>::CancerGenotype(G&& germline, const std::shared_ptr<MappableType>& somatic)
: germline_ {std::forward<G>(germline)}
, somatic_ {somatic}
{}

template <typename MappableType>
template <typename G, typename S>
CancerGenotype<MappableType>::CancerGenotype(G&& germline, S&& somatic)
: germline_ {std::forward<G>(germline)}
, somatic_ {std::forward<S>(somatic)}
{}

template <typename MappableType>
const GenomicRegion& CancerGenotype<MappableType>::mapped_region() const noexcept
{
    using octopus::mapped_region;
    return mapped_region(germline_);
}

template <typename MappableType>
const MappableType& CancerGenotype<MappableType>::operator[](unsigned n) const
{
    return (n < germline_ploidy()) ? germline_[n] : somatic_[n];
}

template <typename MappableType>
const Genotype<MappableType>& CancerGenotype<MappableType>::germline() const
{
    return germline_;
}

template <typename MappableType>
const Genotype<MappableType>& CancerGenotype<MappableType>::somatic() const
{
    return somatic_;
}

template <typename MappableType>
unsigned CancerGenotype<MappableType>::germline_ploidy() const noexcept
{
    return germline_.ploidy();
}

template <typename MappableType>
unsigned CancerGenotype<MappableType>::somatic_ploidy() const noexcept
{
    return somatic_.ploidy();
}

template <typename MappableType>
unsigned CancerGenotype<MappableType>::ploidy() const noexcept
{
    return germline_ploidy() + somatic_ploidy();
}

template <typename MappableType>
bool CancerGenotype<MappableType>::contains(const MappableType& element) const
{
    return germline_.contains(element) || somatic_.contains(element);
}

template <typename MappableType>
unsigned CancerGenotype<MappableType>::count(const MappableType& element) const
{
    return germline_.count(element) + somatic_.count(element);
}

template <typename MappableType>
bool CancerGenotype<MappableType>::is_homozygous() const
{
    return germline_.is_homozygous() && somatic_.count(germline_[0]) == somatic_.ploidy();
}

template <typename MappableType>
unsigned CancerGenotype<MappableType>::zygosity() const
{
    if (somatic_.ploidy() == 1) {
        return germline_.zygosity() + ((germline_.contains(somatic_)) ? 0 : 1);
    } else {
        return copy_unique().size();
    }
}

template <typename MappableType>
std::vector<MappableType> CancerGenotype<MappableType>::copy_unique() const
{
    auto result = germline_.copy_unique();
    auto itr = utils::append(somatic_.copy_unique(), result);
    std::inplace_merge(std::begin(result), itr, std::end(result));
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

// non-member methods

bool contains(const CancerGenotype<Haplotype>& genotype, const Allele& allele);
bool includes(const CancerGenotype<Haplotype>& genotype, const Allele& allele);

template <typename MappableType2, typename MappableType1>
CancerGenotype<MappableType2> copy(const CancerGenotype<MappableType1>& genotype, const GenomicRegion& region)
{
    return CancerGenotype<MappableType2> {copy<MappableType2>(genotype.germline(), region),
                                          copy<MappableType2>(genotype.somatic(), region)};
}

template <typename MappableType1, typename MappableType2>
bool contains(const CancerGenotype<MappableType1>& lhs, const CancerGenotype<MappableType2>& rhs)
{
    return copy<MappableType2>(lhs, rhs.mapped_region()) == rhs;
}

struct CancerGenotypeIndex
{
    GenotypeIndex germline, somatic;
};

std::vector<CancerGenotype<Haplotype>>
generate_all_cancer_genotypes(const std::vector<Genotype<Haplotype>>& germline_genotypes,
                              const std::vector<Haplotype>& somatic_haplotypes,
                              unsigned somatic_ploidy = 1, bool allow_shared = false);

std::vector<CancerGenotype<Haplotype>>
generate_all_cancer_genotypes(const std::vector<Genotype<Haplotype>>& germline_genotypes,
                              const std::vector<GenotypeIndex>& germline_genotype_indices,
                              const std::vector<Haplotype>& somatic_haplotypes,
                              std::vector<CancerGenotypeIndex>& cancer_genotype_indices,
                              unsigned somatic_ploidy = 1, bool allow_shared = false);

template <typename MappableType>
Genotype<MappableType> demote(const CancerGenotype<MappableType>& genotype)
{
    Genotype<MappableType> result {genotype.ploidy()};
    for (const auto& e : genotype.germline()) {
        result.emplace(e);
    }
    for (const auto& e : genotype.somatic()) {
        result.emplace(e);
    }
    return result;
}

template <typename MappableType>
bool operator==(const CancerGenotype<MappableType>& lhs, const CancerGenotype<MappableType>& rhs)
{
    return lhs.germline() == rhs.germline() && lhs.somatic() == rhs.somatic();
}

template <typename MappableType>
std::ostream& operator<<(std::ostream& os, const CancerGenotype<MappableType>& genotype)
{
    os << genotype.germline() << "," << genotype.somatic() << "(cancer)";
    return os;
}

struct CancerGenotypeHash
{
    template <typename T>
    std::size_t operator()(const CancerGenotype<T>& genotype) const
    {
        using boost::hash_combine;
        size_t result {};
        hash_combine(result, std::hash<Genotype<T>>()(genotype.germline()));
        hash_combine(result, std::hash<Genotype<T>>()(genotype.somatic()));
        return result;
    }
};

namespace debug {

template <typename S>
void print_alleles(S&& stream, const CancerGenotype<Haplotype>& genotype)
{
    print_alleles(stream, genotype.germline());
    stream << " + ";
    print_alleles(stream, genotype.somatic());
}

void print_alleles(const CancerGenotype<Haplotype>& genotype);

template <typename S>
void print_variant_alleles(S&& stream, const CancerGenotype<Haplotype>& genotype)
{
    print_variant_alleles(stream, genotype.germline());
    stream << " + ";
    print_variant_alleles(stream, genotype.somatic());
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
