//
//  variant.h
//  Octopus
//
//  Created by Daniel Cooke on 05/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_variant_h
#define Octopus_variant_h

#include <string>
#include <cstdint>
#include <cstddef>
#include <functional> // std::function

#include "genomic_region.h"
#include "comparable.h"

using std::uint_fast32_t;
using std::size_t;

/**
   All variants are considered replacements, i.e. a region of the reference contig is removed,
   and replaced with a novel sequence.
 
   The underlying logic is designed to be transparent to the variant type. The exception is
   modelling the prior probability of different variants. This solution uses a strategy pattern
   injection (the functional 'prior_model' is injected) to acheive runtime polymorphism.
 */
class Variant : Comparable<Variant>
{
public:
    using SizeType   = GenomicRegion::SizeType;
    using StringType = std::string;
    
    Variant() = delete;
    template <typename T1, typename T2>
    explicit Variant(GenomicRegion the_reference_allele_region, T1&& the_reference_allele,
                     T2&& the_alternative_allele, std::function<double()> prior_model);
    
    Variant(const Variant&)            = default;
    Variant& operator=(const Variant&) = default;
    Variant(Variant&&)                 = default;
    Variant& operator=(Variant&&)      = default;
    
    const GenomicRegion& get_reference_allele_region() const noexcept;
    const StringType& get_reference_allele() const noexcept;
    const StringType& get_alternative_allele() const noexcept;
    double get_prior_probability() const noexcept;
    
private:
    GenomicRegion the_reference_allele_region_;
    StringType the_reference_allele;
    StringType the_alternative_allele;
    std::function<double()> prior_model_;
};

template <typename T1, typename T2>
Variant::Variant(GenomicRegion the_reference_allele_region, T1&& the_reference_allele,
                 T2&& the_alternative_allele, std::function<double()> prior_model)
:the_reference_allele_region_ {std::move(the_reference_allele_region)},
the_reference_allele {std::forward<T1>(the_reference_allele)},
the_alternative_allele {std::forward<T2>(the_alternative_allele)},
prior_model_ {prior_model}
{}

inline const GenomicRegion& Variant::get_reference_allele_region() const noexcept
{
    return the_reference_allele_region_;
}

inline const Variant::StringType& Variant::get_reference_allele() const noexcept
{
    return the_reference_allele;
}

inline const Variant::StringType& Variant::get_alternative_allele() const noexcept
{
    return the_alternative_allele;
}

inline double Variant::get_prior_probability() const noexcept
{
    return prior_model_();
}

inline bool overlaps(const Variant& lhs, const Variant& rhs) noexcept
{
    return overlaps(lhs.get_reference_allele_region(), rhs.get_reference_allele_region());
}

inline Variant::SizeType reference_allele_size(const Variant& a_variant) noexcept
{
    return size(a_variant.get_reference_allele_region());
}

inline Variant::SizeType alternative_allele_size(const Variant& a_variant) noexcept
{
    return static_cast<Variant::SizeType>(a_variant.get_alternative_allele().size());
}

inline GenomicRegion::StringType reference_contig_name(const Variant& a_variant)
{
    return a_variant.get_reference_allele_region().get_contig_name();
}

inline GenomicRegion::SizeType reference_allele_begin(const Variant& a_variant)
{
    return a_variant.get_reference_allele_region().get_begin();
}

inline GenomicRegion::SizeType reference_allele_end(const Variant& a_variant)
{
    return a_variant.get_reference_allele_region().get_end();
}

inline bool operator==(const Variant& lhs, const Variant& rhs)
{
    return lhs.get_reference_allele_region() == rhs.get_reference_allele_region() &&
           lhs.get_reference_allele() == rhs.get_reference_allele() &&
           lhs.get_alternative_allele() == rhs.get_alternative_allele();
}

inline bool operator<(const Variant& lhs, const Variant& rhs)
{
    return lhs.get_reference_allele_region() < rhs.get_reference_allele_region();
}

namespace std {
    template <> struct hash<Variant>
    {
        size_t operator()(const Variant& v) const
        {
            return hash<GenomicRegion>()(v.get_reference_allele_region());
        }
    };
}

#endif
