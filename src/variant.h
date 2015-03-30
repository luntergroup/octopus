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
#include <ostream>

#include "genomic_region.h"
#include "comparable.h"

using std::uint_fast32_t;
using std::size_t;

class Variant : public Comparable<Variant>
{
public:
    using SizeType     = GenomicRegion::SizeType;
    using SequenceType = std::string;
    
    Variant() = default;
    
    template <typename GenomicRegion_, typename StringType1, typename StringType2>
    explicit Variant(GenomicRegion_&& the_reference_allele_region, StringType1&& the_reference_allele,
                     StringType2&& the_alternative_allele, std::function<double()> prior_model);
    
    template <typename StringType1, typename StringType2, typename StringType3>
    explicit Variant(StringType1&& the_reference_contig_name, SizeType the_reference_begin,
                     StringType2&& the_reference_allele, StringType3&& the_alternative_allele,
                     std::function<double()> prior_model);
    
    Variant(const Variant&)            = default;
    Variant& operator=(const Variant&) = default;
    Variant(Variant&&)                 = default;
    Variant& operator=(Variant&&)      = default;
    
    const GenomicRegion& get_reference_allele_region() const noexcept;
    const SequenceType& get_reference_allele() const noexcept;
    const SequenceType& get_alternative_allele() const noexcept;
    double get_prior_probability() const noexcept;
    
private:
    // Don't change the order of these members
    SequenceType the_reference_allele_;
    GenomicRegion the_reference_allele_region_;
    SequenceType the_alternative_allele_;
    std::function<double()> prior_model_;
};

template <typename GenomicRegion_, typename StringType1, typename StringType2>
Variant::Variant(GenomicRegion_&& the_reference_allele_region, StringType1&& the_reference_allele,
                 StringType2&& the_alternative_allele, std::function<double()> prior_model)
:
the_reference_allele_ {std::forward<StringType1>(the_reference_allele)},
the_reference_allele_region_ {std::forward<GenomicRegion_>(the_reference_allele_region)},
the_alternative_allele_ {std::forward<StringType2>(the_alternative_allele)},
prior_model_ {prior_model}
{}

template <typename StringType1, typename StringType2, typename StringType3>
Variant::Variant(StringType1&& the_reference_contig_name, SizeType the_reference_begin,
                 StringType2&& the_reference_allele, StringType3&& the_alternative_allele,
                 std::function<double()> prior_model)
:
// The reference allele needs to be initialised first so it's size can be used for the
// region initialisation.
the_reference_allele_ {std::forward<StringType2>(the_reference_allele)},
the_reference_allele_region_ {GenomicRegion{std::forward<StringType1>(the_reference_contig_name),
    the_reference_begin, the_reference_begin + static_cast<SizeType>(the_reference_allele_.size())}},
the_alternative_allele_ {std::forward<StringType3>(the_alternative_allele)},
prior_model_ {prior_model}
{}

inline const GenomicRegion& Variant::get_reference_allele_region() const noexcept
{
    return the_reference_allele_region_;
}

inline const Variant::SequenceType& Variant::get_reference_allele() const noexcept
{
    return the_reference_allele_;
}

inline const Variant::SequenceType& Variant::get_alternative_allele() const noexcept
{
    return the_alternative_allele_;
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

inline std::ostream& operator<<(std::ostream& os, const Variant& a_variant)
{
    os << a_variant.get_reference_allele_region() << " " << a_variant.get_reference_allele() << " " <<
        a_variant.get_alternative_allele();
    return os;
}

#endif
