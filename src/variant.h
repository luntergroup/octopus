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
#include <ostream>

#include "genomic_region.h"
#include "comparable.h"
#include "mappable.h"

using std::uint_fast32_t;

class Variant : public Comparable<Variant>, public Mappable<Variant>
{
public:
    using SizeType     = GenomicRegion::SizeType;
    using SequenceType = std::string;
    
    Variant() = default;
    
    template <typename GenomicRegion_, typename StringType1, typename StringType2>
    explicit Variant(GenomicRegion_&& the_reference_allele_region, StringType1&& the_reference_allele,
                     StringType2&& the_alternative_allele, double the_segregation_probability);
    
    template <typename StringType1, typename StringType2, typename StringType3>
    explicit Variant(StringType1&& the_reference_contig_name, SizeType the_reference_begin,
                     StringType2&& the_reference_allele, StringType3&& the_alternative_allele,
                     double the_segregation_probability);
    ~Variant() = default;
    
    Variant(const Variant&)            = default;
    Variant& operator=(const Variant&) = default;
    Variant(Variant&&)                 = default;
    Variant& operator=(Variant&&)      = default;
    
    const GenomicRegion& get_region() const noexcept;
    const SequenceType& get_reference_allele() const noexcept;
    const SequenceType& get_alternative_allele() const noexcept;
    void set_segregation_probability(double p) noexcept;
    double get_segregation_probability() const noexcept;
    
private:
    //WARNING: Don't change the order of these members! Order required by second constructor.
    SequenceType the_reference_allele_;
    GenomicRegion the_reference_allele_region_;
    SequenceType the_alternative_allele_;
    double the_segregation_probability_;
};

template <typename GenomicRegion_, typename StringType1, typename StringType2>
Variant::Variant(GenomicRegion_&& the_reference_allele_region, StringType1&& the_reference_allele,
                 StringType2&& the_alternative_allele, double the_segregation_probability)
:
the_reference_allele_ {std::forward<StringType1>(the_reference_allele)},
the_reference_allele_region_ {std::forward<GenomicRegion_>(the_reference_allele_region)},
the_alternative_allele_ {std::forward<StringType2>(the_alternative_allele)},
the_segregation_probability_ {the_segregation_probability}
{}

template <typename StringType1, typename StringType2, typename StringType3>
Variant::Variant(StringType1&& the_reference_contig_name, SizeType the_reference_begin,
                 StringType2&& the_reference_allele, StringType3&& the_alternative_allele,
                 double the_segregation_probability)
:
// The reference allele needs to be initialised first so it's size can be used for the
// region initialisation.
the_reference_allele_ {std::forward<StringType2>(the_reference_allele)},
the_reference_allele_region_ {GenomicRegion{std::forward<StringType1>(the_reference_contig_name),
    the_reference_begin, the_reference_begin + static_cast<SizeType>(the_reference_allele_.size())}},
the_alternative_allele_ {std::forward<StringType3>(the_alternative_allele)},
the_segregation_probability_ {the_segregation_probability}
{}

inline const GenomicRegion& Variant::get_region() const noexcept
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

inline void Variant::set_segregation_probability(double p) noexcept
{
    the_segregation_probability_ = p;
}

inline double Variant::get_segregation_probability() const noexcept
{
    return the_segregation_probability_;
}

inline Variant::SizeType reference_allele_size(const Variant& a_variant) noexcept
{
    return size(a_variant.get_region());
}

inline Variant::SizeType alternative_allele_size(const Variant& a_variant) noexcept
{
    return static_cast<Variant::SizeType>(a_variant.get_alternative_allele().size());
}

inline GenomicRegion::StringType reference_contig_name(const Variant& a_variant)
{
    return a_variant.get_region().get_contig_name();
}

inline GenomicRegion::SizeType reference_allele_begin(const Variant& a_variant)
{
    return a_variant.get_region().get_begin();
}

inline GenomicRegion::SizeType reference_allele_end(const Variant& a_variant)
{
    return a_variant.get_region().get_end();
}

inline bool operator==(const Variant& lhs, const Variant& rhs)
{
    return lhs.get_region() == rhs.get_region() &&
           lhs.get_reference_allele() == rhs.get_reference_allele() &&
           lhs.get_alternative_allele() == rhs.get_alternative_allele();
}

inline bool operator<(const Variant& lhs, const Variant& rhs)
{
    // This check is required for consistency with operator==
    if (lhs.get_region() == rhs.get_region()) {
        return (lhs.get_reference_allele() < rhs.get_reference_allele()) ? true :
                            lhs.get_alternative_allele() < rhs.get_alternative_allele();
    } else {
        return lhs.get_region() < rhs.get_region();
    }
}

namespace std {
    template <> struct hash<Variant>
    {
        size_t operator()(const Variant& v) const
        {
            return hash<GenomicRegion>()(v.get_region());
        }
    };
}

inline std::ostream& operator<<(std::ostream& os, const Variant& a_variant)
{
    os << a_variant.get_region() << " " << a_variant.get_reference_allele() << " " <<
        a_variant.get_alternative_allele();
    return os;
}

#endif
