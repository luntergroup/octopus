//
//  variant.hpp
//  Octopus
//
//  Created by Daniel Cooke on 05/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_variant_hpp
#define Octopus_variant_hpp

#include <string>
#include <ostream>

#include <boost/functional/hash.hpp>

#include "genomic_region.hpp"
#include "allele.hpp"
#include "comparable.hpp"
#include "mappable.hpp"

/*
 A variant is a combination of a reference allele and an alternative allele.
 */
class Variant : public Comparable<Variant>, public Mappable<Variant>
{
public:
    using SizeType     = Allele::SizeType;
    using SequenceType = Allele::SequenceType;
    
    Variant() = default;
    
    explicit Variant(Allele reference, Allele alternative);
    
    template <typename GenomicRegion_, typename SequenceType1, typename SequenceType2>
    explicit Variant(GenomicRegion_&& reference_allele_region, SequenceType1&& reference_allele,
                     SequenceType2&& alternative_allele);
    
    template <typename SequenceType1, typename SequenceType2, typename SequenceType3>
    explicit Variant(SequenceType1&& reference_contig_name, SizeType reference_begin,
                     SequenceType2&& reference_allele, SequenceType3&& alternative_allele);
    
    ~Variant() = default;
    
    Variant(const Variant&)            = default;
    Variant& operator=(const Variant&) = default;
    Variant(Variant&&)                 = default;
    Variant& operator=(Variant&&)      = default;
    
    const GenomicRegion& get_region() const noexcept;
    
    const Allele& get_ref_allele() const noexcept;
    const Allele& get_alt_allele() const noexcept;
    
private:
    Allele ref_, alt_;
};

template <typename GenomicRegion_, typename SequenceType1, typename SequenceType2>
Variant::Variant(GenomicRegion_&& ref_region, SequenceType1&& ref_sequence,
                 SequenceType2&& alt_sequence)
:
ref_ {ref_region, std::forward<SequenceType1>(ref_sequence)},
alt_ {std::forward<GenomicRegion_>(ref_region), std::forward<SequenceType2>(alt_sequence)}
{}

template <typename SequenceType1, typename SequenceType2, typename SequenceType3>
Variant::Variant(SequenceType1&& ref_contig_name, SizeType ref_begin,
                 SequenceType2&& ref_sequence, SequenceType3&& alt_sequence)
:
ref_ {std::forward<SequenceType1>(ref_contig_name), ref_begin, std::forward<SequenceType2>(ref_sequence)},
alt_ {ref_.get_region(), std::forward<SequenceType3>(alt_sequence)}
{}

// non-member methods

const Variant::SequenceType& get_ref_sequence(const Variant& variant);
const Variant::SequenceType& get_alt_sequence(const Variant& variant);

Variant::SizeType ref_sequence_size(const Variant& variant);
Variant::SizeType alt_sequence_size(const Variant& variant);

bool operator==(const Variant& lhs, const Variant& rhs);
bool operator<(const Variant& lhs, const Variant& rhs);

namespace std {
    template <> struct hash<Variant>
    {
        size_t operator()(const Variant& variant) const
        {
            using boost::hash_combine;
            size_t result {};
            hash_combine(result, hash<GenomicRegion>()(variant.get_region()));
            hash_combine(result, hash<Allele::SequenceType>()(get_ref_sequence(variant)));
            hash_combine(result, hash<Allele::SequenceType>()(get_alt_sequence(variant)));
            return result;
        }
    };
} // end namespace std

std::ostream& operator<<(std::ostream& os, const Variant& variant);

#endif
