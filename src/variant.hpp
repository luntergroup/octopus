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

class ReferenceGenome;

/*
 A variant is a combination of a reference allele and an alternative allele.
 */
class Variant : public Comparable<Variant>, public Mappable<Variant>
{
public:
    using SizeType     = Allele::SizeType;
    using SequenceType = Allele::SequenceType;
    
    Variant() = default;
    
    explicit Variant(const Allele& reference, const Allele& alternative);
    explicit Variant(Allele&& reference, Allele&& alternative);
    
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
    Allele reference_, alternative_;
};

template <typename GenomicRegion_, typename SequenceType1, typename SequenceType2>
Variant::Variant(GenomicRegion_&& ref_region, SequenceType1&& ref_sequence,
                 SequenceType2&& alt_sequence)
:
reference_ {ref_region, std::forward<SequenceType1>(ref_sequence)},
alternative_ {std::forward<GenomicRegion_>(ref_region), std::forward<SequenceType2>(alt_sequence)}
{}

template <typename SequenceType1, typename SequenceType2, typename SequenceType3>
Variant::Variant(SequenceType1&& ref_contig_name, SizeType ref_begin,
                 SequenceType2&& ref_sequence, SequenceType3&& alt_sequence)
:
reference_ {std::forward<SequenceType1>(ref_contig_name), ref_begin, std::forward<SequenceType2>(ref_sequence)},
alternative_ {reference_.get_region(), std::forward<SequenceType3>(alt_sequence)}
{}

// non-member methods

Variant make_variant(const Allele& alt_allele, const ReferenceGenome& reference);
Variant make_variant(const std::string& region_str, Variant::SequenceType alt_sequence,
                     const ReferenceGenome& reference);

const Variant::SequenceType& get_ref_sequence(const Variant& variant);
const Variant::SequenceType& get_alt_sequence(const Variant& variant);

Variant::SizeType ref_sequence_size(const Variant& variant);
Variant::SizeType alt_sequence_size(const Variant& variant);

bool operator==(const Variant& lhs, const Variant& rhs);
bool operator<(const Variant& lhs, const Variant& rhs);

void remove_duplicates(std::vector<Variant>& variants);

/*
 Decomposes a list of variants into unique alleles in the same order as the given variants.
 The reference allele for each unique GenomicRegion will appear first in the sub-list.
 */
std::vector<Allele> decompose(const std::vector<Variant>& variants);

std::vector<Allele> get_intervening_reference_alleles(const std::vector<Variant>& variants,
                                                      const ReferenceGenome& reference);

/*
 A variant is parsimonious if and only if it is represented in as few nucleotides as possible
 without an allele of length 0.
 */
bool is_parsimonious(const Variant& variant) noexcept;

Variant make_parsimonious(const Variant& variant, const ReferenceGenome& reference);

bool is_left_alignable(const Variant& variant) noexcept;

// Note there are no 'is_left_aligned' or 'is_normalised' methods because they
// would take as much work as calling 'left_align'.

/*
 A variant is left aligned if and only if it is no longer possible to shift its position
 to the left while keeping the length of all its alleles constant.
 */
Variant left_align(const Variant& variant, const ReferenceGenome& reference,
                   Variant::SizeType extension_size = 30);

/*
 A variant is normalised if and only if it is parsimonious and left aligned.
 */
Variant normalise(const Variant& variant, const ReferenceGenome& reference,
                  Variant::SizeType extension_size = 30);

Variant pad_left(const Variant& variant, const Variant::SequenceType& sequence);
Variant pad_right(const Variant& variant, const Variant::SequenceType& sequence);
Variant pad_right(const Variant& variant, const ReferenceGenome& reference, Variant::SizeType n);
Variant pad_left(const Variant& variant, const ReferenceGenome& reference, Variant::SizeType n);
Variant pad_right(const Variant& variant, const ReferenceGenome& reference, Variant::SizeType n);

/*
 Left aligns all input Variants and removes any resulting duplicates. The returned variants are sorted.
 */
std::vector<Variant> unique_left_align(const std::vector<Variant>& variants, const ReferenceGenome& reference);

std::vector<Variant> parsimonise_each(const std::vector<Variant>& variants, const ReferenceGenome& reference);

std::vector<Variant> parsimonise_together(const std::vector<Variant>& variants, const ReferenceGenome& reference);

bool is_snp(const Variant& variant) noexcept;

bool is_insertion(const Variant& variant) noexcept;

bool is_deletion(const Variant& variant) noexcept;

bool is_indel(const Variant& variant) noexcept;

bool is_mnv(const Variant& variant) noexcept;

bool is_transition(const Variant& variant) noexcept;

bool is_transversion(const Variant& variant) noexcept;

std::vector<Allele::SequenceType> get_alt_allele_sequences(const std::vector<Variant>& variants);

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
