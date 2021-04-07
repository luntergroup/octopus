// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef cancer_genotype_hpp
#define cancer_genotype_hpp

#include <initializer_list>
#include <utility>
#include <memory>
#include <functional>
#include <iterator>
#include <algorithm>
#include <cassert>
#include <iostream>

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
    
    CancerGenotype(const CancerGenotype&)            = default;
    CancerGenotype& operator=(const CancerGenotype&) = default;
    CancerGenotype(CancerGenotype&&)                 = default;
    CancerGenotype& operator=(CancerGenotype&&)      = default;
    
    ~CancerGenotype() = default;
    
    const GenomicRegion& mapped_region() const noexcept;
    
    const MappableType& operator[](unsigned n) const noexcept;
    MappableType& operator[](unsigned n) noexcept;
    
    const Genotype<MappableType>& germline() const noexcept;
    Genotype<MappableType>& germline() noexcept;
    const Genotype<MappableType>& somatic() const noexcept;
    Genotype<MappableType>& somatic() noexcept;
    
    unsigned germline_ploidy() const noexcept;
    unsigned somatic_ploidy() const noexcept;
    unsigned ploidy() const noexcept;
    
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
const MappableType& CancerGenotype<MappableType>::operator[](unsigned n) const noexcept
{
    return (n < germline_ploidy()) ? germline_[n] : somatic_[n];
}

template <typename MappableType>
MappableType& CancerGenotype<MappableType>::operator[](unsigned n) noexcept
{
    return (n < germline_ploidy()) ? germline_[n] : somatic_[n];
}

template <typename MappableType>
const Genotype<MappableType>& CancerGenotype<MappableType>::germline() const noexcept
{
    return germline_;
}

template <typename MappableType>
Genotype<MappableType>& CancerGenotype<MappableType>::germline() noexcept
{
    return germline_;
}

template <typename MappableType>
const Genotype<MappableType>& CancerGenotype<MappableType>::somatic() const noexcept
{
    return somatic_;
}

template <typename MappableType>
Genotype<MappableType>& CancerGenotype<MappableType>::somatic() noexcept
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

// free functions

template <typename MappableType>
auto count(const CancerGenotype<MappableType>& genotype)
{
    return count(genotype.germline()) + genotype(genotype.somatic());
}

template <typename MappableType>
bool is_homozygous(const CancerGenotype<MappableType>& genotype)
{
    return is_homozygous(genotype.germline()) && count(genotype.somatic(), genotype.germline()[0]) == genotype.somatic_ploidy();
}

template <typename MappableType>
auto zygosity(const CancerGenotype<MappableType>& genotype)
{
    if (genotype.somatic_ploidy() == 1) {
        return zygosity(genotype.germline()) + ((contains(genotype.germline(), genotype.somatic()[0])) ? 0 : 1);
    } else {
        auto result = zygosity(genotype.germline());
        for (const auto& element : collapse(genotype.somatic())) {
            if (!contains(genotype.germline()), element) {
                ++result;
            }
        }
        return result;
    }
}

template <typename MappableType>
CancerGenotype<MappableType> collapse(const CancerGenotype<MappableType>& genotype)
{
    return {collapse(genotype.germline()), collapse(genotype.somatic())};
}

// non-member methods

template <typename MappableType1, typename MappableType2>
bool contains(const CancerGenotype<MappableType1>& genotype, const MappableType2& element)
{
    return contains(genotype.germline(), element) || contains(genotype.somatic(), element);
}

template <typename MappableType1, typename MappableType2>
bool includes(const CancerGenotype<MappableType1>& genotype, const MappableType2& element)
{
    return includes(genotype.germline(), element) || includes(genotype.somatic(), element);
}

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

template <typename MappableType>
MappableBlock<CancerGenotype<MappableType>>
generate_all_cancer_genotypes(const MappableBlock<Genotype<MappableType>>& germline_genotypes,
                              const MappableBlock<MappableType>& elements,
                              unsigned somatic_ploidy = 1,
                              bool allow_shared = false)
{
    const auto somatic_genotypes = generate_all_max_zygosity_genotypes(elements, somatic_ploidy);
    MappableBlock<CancerGenotype<MappableType>> result {mapped_region(elements)};
    result.reserve(germline_genotypes.size() * somatic_genotypes.size());
    for (const auto& germline : germline_genotypes) {
        for (const auto& somatic : somatic_genotypes) {
            if (allow_shared || !have_shared(germline, somatic)) {
                result.emplace_back(germline, somatic);
            }
        }
    }
    return result;
}

template <typename MappableType>
MappableBlock<CancerGenotype<MappableType>>
extend_somatic(const MappableBlock<CancerGenotype<MappableType>>& old_genotypes,
               const MappableBlock<MappableType>& elements,
               bool allow_shared = false)
{
    MappableBlock<CancerGenotype<MappableType>> result {mapped_region(elements)};
    result.reserve(old_genotypes.size() * elements.size());
    for (const auto& old_genotype : old_genotypes) {
        for (const auto& element : elements) {
            if (allow_shared || !contains(old_genotype.germline(), element)) {
                auto new_somatic_genotype = old_genotype.somatic();
                new_somatic_genotype.emplace(element);
                if (is_max_zygosity(new_somatic_genotype)) {
                    result.emplace_back(old_genotype.germline(), std::move(new_somatic_genotype));
                }
            }
        }
    }
    return result;
}

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

template <typename S, typename T>
void print_alleles(S&& stream, const CancerGenotype<T>& genotype)
{
    print_alleles(stream, genotype.germline());
    stream << " + ";
    print_alleles(stream, genotype.somatic());
}

template <typename S, typename T>
void print_variant_alleles(S&& stream, const CancerGenotype<T>& genotype)
{
    print_variant_alleles(stream, genotype.germline());
    stream << " + ";
    print_variant_alleles(stream, genotype.somatic());
}

template <typename T>
void print_alleles(const CancerGenotype<T>& genotype)
{
    print_alleles(std::cout, genotype);
}

template <typename T>
void print_variant_alleles(const CancerGenotype<T>& genotype)
{
    print_variant_alleles(std::cout, genotype);
}

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
