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
#include <ostream>

#include "genomic_region.h"
#include "allele.h"
#include "comparable.h"
#include "mappable.h"

/**
 A variant is a combination of a reference allele and an alternative allele.
 */
class Variant : public Comparable<Variant>, public Mappable<Variant>
{
public:
    using SizeType     = GenomicRegion::SizeType;
    using SequenceType = Allele::SequenceType;
    using RealType     = Allele::RealType;
    
    Variant() = default;
    template <typename GenomicRegion_, typename SequenceType1, typename SequenceType2>
    explicit Variant(GenomicRegion_&& the_reference_allele_region, SequenceType1&& the_reference_allele,
                     SequenceType2&& the_alternative_allele, RealType the_reference_allele_probability,
                     RealType the_alternative_allele_probability);
    template <typename SequenceType1, typename SequenceType2, typename SequenceType3>
    explicit Variant(SequenceType1&& the_reference_contig_name, SizeType the_reference_begin,
                     SequenceType2&& the_reference_allele, SequenceType3&& the_alternative_allele,
                     RealType the_reference_allele_probability, RealType the_alternative_allele_probability);
    ~Variant() = default;
    
    Variant(const Variant&)            = default;
    Variant& operator=(const Variant&) = default;
    Variant(Variant&&)                 = default;
    Variant& operator=(Variant&&)      = default;
    
    const GenomicRegion& get_region() const noexcept;
    
    Allele get_reference_allele() const;
    Allele get_alternative_allele() const;
    
    const SequenceType& get_reference_allele_sequence() const noexcept;
    const SequenceType& get_alternative_allele_sequence() const noexcept;
    
    void set_reference_allele_probability(RealType probability) noexcept;
    void set_alternative_allele_probability(RealType probability) noexcept;
    RealType get_reference_allele_probability() const noexcept;
    RealType get_alternative_allele_probability() const noexcept;
    
private:
    //WARNING: Do not change the order of these members! Order required by second constructor.
    SequenceType the_reference_allele_sequence_;
    GenomicRegion the_reference_allele_region_;
    SequenceType the_alternative_allele_sequence_;
    RealType the_reference_allele_probability_;
    RealType the_alternative_allele_probability_;
};

template <typename GenomicRegion_, typename SequenceType1, typename SequenceType2>
Variant::Variant(GenomicRegion_&& the_reference_allele_region, SequenceType1&& the_reference_allele,
                 SequenceType2&& the_alternative_allele, RealType the_reference_allele_probability,
                 RealType the_alternative_allele_probability)
:
the_reference_allele_sequence_ {std::forward<SequenceType1>(the_reference_allele)},
the_reference_allele_region_ {std::forward<GenomicRegion_>(the_reference_allele_region)},
the_alternative_allele_sequence_ {std::forward<SequenceType2>(the_alternative_allele)},
the_reference_allele_probability_ {the_reference_allele_probability},
the_alternative_allele_probability_ {the_alternative_allele_probability}
{}

template <typename SequenceType1, typename SequenceType2, typename SequenceType3>
Variant::Variant(SequenceType1&& the_reference_contig_name, SizeType the_reference_begin,
                 SequenceType2&& the_reference_allele, SequenceType3&& the_alternative_allele,
                 RealType the_reference_allele_probability, RealType the_alternative_allele_probability)
:
// The reference allele needs to be initialised first so it's size can be used for the
// region initialisation.
the_reference_allele_sequence_ {std::forward<SequenceType2>(the_reference_allele)},
the_reference_allele_region_ {GenomicRegion{std::forward<SequenceType1>(the_reference_contig_name),
    the_reference_begin, the_reference_begin + static_cast<SizeType>(the_reference_allele_sequence_.size())}},
the_alternative_allele_sequence_ {std::forward<SequenceType3>(the_alternative_allele)},
the_reference_allele_probability_ {the_reference_allele_probability},
the_alternative_allele_probability_ {the_alternative_allele_probability}
{}

inline const GenomicRegion& Variant::get_region() const noexcept
{
    return the_reference_allele_region_;
}

inline Allele Variant::get_reference_allele() const
{
    return Allele {the_reference_allele_region_, the_reference_allele_sequence_, the_reference_allele_probability_};
}

inline Allele Variant::get_alternative_allele() const
{
    return Allele {the_reference_allele_region_, the_alternative_allele_sequence_, the_alternative_allele_probability_};
}

inline const Variant::SequenceType& Variant::get_reference_allele_sequence() const noexcept
{
    return the_reference_allele_sequence_;
}

inline const Variant::SequenceType& Variant::get_alternative_allele_sequence() const noexcept
{
    return the_alternative_allele_sequence_;
}

inline void Variant::set_reference_allele_probability(RealType probability) noexcept
{
    the_reference_allele_probability_ = probability;
}

inline void Variant::set_alternative_allele_probability(RealType probability) noexcept
{
    the_alternative_allele_probability_ = probability;
}

inline Variant::RealType Variant::get_reference_allele_probability() const noexcept
{
    return the_reference_allele_probability_;
}

inline Variant::RealType Variant::get_alternative_allele_probability() const noexcept
{
    return the_alternative_allele_probability_;
}

inline Variant::SizeType reference_allele_size(const Variant& a_variant) noexcept
{
    return size(a_variant.get_region());
}

inline Variant::SizeType alternative_allele_size(const Variant& a_variant) noexcept
{
    return static_cast<Variant::SizeType>(a_variant.get_alternative_allele_sequence().size());
}

inline bool operator==(const Variant& lhs, const Variant& rhs)
{
    return lhs.get_region() == rhs.get_region() &&
           lhs.get_reference_allele_sequence() == rhs.get_reference_allele_sequence() &&
           lhs.get_alternative_allele_sequence() == rhs.get_alternative_allele_sequence();
}

inline bool operator<(const Variant& lhs, const Variant& rhs)
{
    // This check is required for consistency with operator==
    if (lhs.get_region() == rhs.get_region()) {
        return (lhs.get_reference_allele_sequence() < rhs.get_reference_allele_sequence()) ? true :
                    lhs.get_alternative_allele_sequence() < rhs.get_alternative_allele_sequence();
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
    os << a_variant.get_region() << " " << a_variant.get_reference_allele_sequence() << " " <<
        a_variant.get_alternative_allele_sequence();
    return os;
}

#endif
