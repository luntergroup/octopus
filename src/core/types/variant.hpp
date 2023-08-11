// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef variant_hpp
#define variant_hpp

#include <string>
#include <functional>
#include <algorithm>
#include <iterator>
#include <utility>
#include <iosfwd>
#include <stdexcept>

#include <boost/functional/hash.hpp>

#include "basics/genomic_region.hpp"
#include "concepts/comparable.hpp"
#include "concepts/mappable.hpp"
#include "allele.hpp"

namespace octopus {

class ReferenceGenome;

/*
 A variant is a combination of a reference allele and an alternative allele.
 */
class Variant : public Comparable<Variant>, public Mappable<Variant>
{
public:
    using MappingDomain      = Allele::MappingDomain;
    using NucleotideSequence = Allele::NucleotideSequence;
    
    Variant() = default;
    
    template <typename R, typename A>
    Variant(R&& reference, A&& alternative);
    
    template <typename GenomicRegion_, typename Sequence1, typename Sequence2>
    Variant(GenomicRegion_&& reference_allele_region,
            Sequence1&& reference_allele, Sequence2&& alternative_allele);
    
    template <typename String, typename Sequence1, typename Sequence2>
    Variant(String&& reference_contig_name, MappingDomain::Position reference_begin,
            Sequence1&& reference_allele, Sequence2&& alternative_allele);
    
    Variant(const Variant&)            = default;
    Variant& operator=(const Variant&) = default;
    Variant(Variant&&)                 = default;
    Variant& operator=(Variant&&)      = default;
    
    ~Variant() = default;
    
    const GenomicRegion& mapped_region() const noexcept;
    
    const Allele& ref_allele() const noexcept;
    const Allele& alt_allele() const noexcept;
    
private:
    Allele reference_, alternative_;
};

template <typename R, typename A>
Variant::Variant(R&& reference, A&& alternative)
: reference_ {std::forward<R>(reference)}
, alternative_ {std::forward<A>(alternative)}
{
    if (!is_same_region(reference_, alternative_)) {
        throw std::logic_error {"Variant: reference & alternative alleles must define the same region"};
    }
}

template <typename GenomicRegion_, typename Sequence1, typename Sequence2>
Variant::Variant(GenomicRegion_&& ref_region, Sequence1&& ref_sequence, Sequence2&& alt_sequence)
: reference_ {ref_region, std::forward<Sequence1>(ref_sequence)}
, alternative_ {std::forward<GenomicRegion_>(ref_region), std::forward<Sequence2>(alt_sequence)}
{}

template <typename String, typename Sequence1, typename Sequence2>
Variant::Variant(String&& ref_contig_name, const MappingDomain::Position ref_begin,
                 Sequence1&& ref_sequence, Sequence2&& alt_sequence)
: reference_ {std::forward<String>(ref_contig_name), ref_begin, std::forward<Sequence1>(ref_sequence)}
, alternative_ {reference_.mapped_region(), std::forward<Sequence2>(alt_sequence)}
{}

// non-member methods

Variant make_variant(const Allele& alt_allele, const ReferenceGenome& reference);

const Variant::NucleotideSequence& ref_sequence(const Variant& variant);
const Variant::NucleotideSequence& alt_sequence(const Variant& variant);

Variant::NucleotideSequence::size_type ref_sequence_size(const Variant& variant);
Variant::NucleotideSequence::size_type alt_sequence_size(const Variant& variant);

bool operator==(const Variant& lhs, const Variant& rhs);
bool operator<(const Variant& lhs, const Variant& rhs);

void remove_duplicates(std::vector<Variant>& variants);

std::vector<Variant> split_mnv(const Variant& variant);

/*
 Decomposes a list of variants into unique alleles in the same order as the given variants.
 The reference allele for each unique GenomicRegion will appear first in the sub-list.
 */
std::vector<Allele> decompose(const std::vector<Variant>& variants);
std::vector<std::reference_wrapper<const Allele>> decompose_ref(const std::vector<Variant>& variants);

std::vector<Allele> extract_intervening_reference_alleles(const std::vector<Variant>& variants,
                                                          const ReferenceGenome& reference);

bool can_trim(const Variant& v);

Variant trim(const Variant& v);

/*
 A variant is parsimonious if and only if it is represented in as few nucleotides as possible
 without an allele of length 0.
 */
bool is_parsimonious(const Variant& variant) noexcept;

/*
 See http://genome.sph.umich.edu/wiki/Variant_Normalization for details on variant parsimony.
 The requirement of the template type R is a SequenceType sequence(const GenomicRegion&)
 method.
 
 G is a functor which must return a char of the given GenomicRegion position.
 */
template <typename G>
Variant make_parsimonious(const Variant& variant, G generator)
{
    using std::cbegin; using std::cend; using std::crbegin; using std::crend;
    
    using NucleotideSequence = Variant::NucleotideSequence;
    
    const auto& old_ref_region   = mapped_region(variant);
    const auto& old_ref_sequence = ref_sequence(variant);
    const auto& old_alt_sequence = alt_sequence(variant);
    
    if (old_ref_sequence.empty() || old_alt_sequence.empty()) {
        if (old_ref_region.begin() > 0) {
            auto new_ref_region = expand_lhs(old_ref_region, 1);
            
            const char new_base = generator(head_position(new_ref_region));
            
            NucleotideSequence new_ref_allele(old_ref_sequence.size() + 1, new_base);
            std::copy(cbegin(old_ref_sequence), cend(old_ref_sequence),
                      std::next(begin(new_ref_allele)));
            
            NucleotideSequence new_alt_allele(old_alt_sequence.size() + 1, new_base);
            std::copy(cbegin(old_alt_sequence), cend(old_alt_sequence),
                      std::next(begin(new_alt_allele)));
            
            return Variant {std::move(new_ref_region), std::move(new_ref_allele), std::move(new_alt_allele)};
        } else {
            // In this rare case the only option is to pad to the right.
            auto new_ref_region = expand_rhs(old_ref_region, 1);
            
            const char new_base = generator(tail_position(new_ref_region));
            
            NucleotideSequence new_ref_allele(old_ref_sequence.size() + 1, new_base);
            std::copy(cbegin(old_ref_sequence), cend(old_ref_sequence), begin(new_ref_allele));
            
            NucleotideSequence new_alt_allele(old_alt_sequence.size() + 1, new_base);
            std::copy(cbegin(old_alt_sequence), cend(old_alt_sequence), begin(new_alt_allele));
            
            return Variant {std::move(new_ref_region), std::move(new_ref_allele), std::move(new_alt_allele)};
        }
    }
    
    const auto& alleles = std::minmax(old_ref_sequence, old_alt_sequence,
                                      [] (const auto& lhs, const auto& rhs) {
                                          return lhs.size() < rhs.size();
                                      });
    const auto& small_allele = alleles.first;
    const auto& big_allele   = alleles.second;
    
    auto pf = std::mismatch(cbegin(small_allele), cend(small_allele), cbegin(big_allele));
    
    if (pf.first == cbegin(small_allele) && small_allele.size() == 1) {
        return variant;
    }
    
    const auto pb = std::mismatch(crbegin(small_allele), std::make_reverse_iterator(pf.first),
                                  crbegin(big_allele));
    
    if (pf.first == cbegin(small_allele) && pb.first == crbegin(small_allele)) {
        return variant;
    }
    
    if (pf.first == std::cend(small_allele)) {
        --pf.first;
        --pf.second;
    }
    
    NucleotideSequence new_big_allele   {pf.second, pb.second.base()};
    NucleotideSequence new_small_allele {pf.first, pb.first.base()};
    
    using S = GenomicRegion::Position;
    GenomicRegion new_ref_region {
        old_ref_region.contig_name(),
        old_ref_region.begin() + static_cast<S>(std::distance(cbegin(small_allele), pf.first)),
        old_ref_region.end()   - static_cast<S>(std::distance(crbegin(small_allele), pb.first))
    };
    
    if (old_ref_sequence.size() > old_alt_sequence.size()) {
        return Variant {std::move(new_ref_region), std::move(new_big_allele), std::move(new_small_allele)};
    } else {
        return Variant {std::move(new_ref_region), std::move(new_small_allele), std::move(new_big_allele)};
    }
}

Variant make_parsimonious(const Variant& variant, const ReferenceGenome& reference);

bool is_left_alignable(const Variant& variant) noexcept;

// Note there are no 'is_left_aligned' or 'is_normalised' methods because they
// would take as much work as calling 'left_align'.

/*
 A variant is left aligned if and only if it is no longer possible to shift its position
 to the left while keeping the length of all its alleles constant.
 */
Variant left_align(const Variant& variant, const ReferenceGenome& reference,
                   GenomicRegion::Size extension_size = 30);

/*
 A variant is normalised if and only if it is parsimonious and left aligned.
 */
Variant normalise(const Variant& variant, const ReferenceGenome& reference,
                  unsigned extension_size = 30);

Variant pad_left(const Variant& variant, const Variant::NucleotideSequence& sequence);
Variant pad_right(const Variant& variant, const Variant::NucleotideSequence& sequence);
Variant pad_right(const Variant& variant, const ReferenceGenome& reference, unsigned n);
Variant pad_left(const Variant& variant, const ReferenceGenome& reference, unsigned n);
Variant pad_right(const Variant& variant, const ReferenceGenome& reference, unsigned n);

/*
 Left aligns all input Variants and removes any resulting duplicates. The returned variants are sorted.
 */
std::vector<Variant> unique_left_align(const std::vector<Variant>& variants,
                                       const ReferenceGenome& reference);
std::vector<Variant> unique_left_align(std::vector<Variant>&& variants,
                                       const ReferenceGenome& reference);

std::vector<Variant> parsimonise_each(const std::vector<Variant>& variants,
                                      const ReferenceGenome& reference);
std::vector<Variant> parsimonise_each(std::vector<Variant>&& variants,
                                      const ReferenceGenome& reference);

std::vector<Variant> parsimonise_together(const std::vector<Variant>& variants,
                                          const ReferenceGenome& reference);

bool is_snv(const Variant& variant) noexcept;
bool is_mnv(const Variant& variant) noexcept;
bool is_insertion(const Variant& variant) noexcept;
bool is_deletion(const Variant& variant) noexcept;
bool is_indel(const Variant& variant) noexcept;
bool is_simple_insertion(const Variant& variant) noexcept;
bool is_simple_deletion(const Variant& variant) noexcept;
bool is_simple_indel(const Variant& variant) noexcept;
bool are_same_type(const Variant& lhs, const Variant& rhs) noexcept;
bool is_transition(const Variant& variant) noexcept;
bool is_transversion(const Variant& variant) noexcept;
Variant::NucleotideSequence::size_type indel_size(const Variant& variant) noexcept;

std::vector<Allele::NucleotideSequence> extract_alt_allele_sequences(const std::vector<Variant>& variants);

std::ostream& operator<<(std::ostream& os, const Variant& variant);

struct VariantHash
{
    std::size_t operator()(const Variant& variant) const
    {
        using boost::hash_combine;
        size_t result {};
        hash_combine(result, std::hash<GenomicRegion>()(variant.mapped_region()));
        hash_combine(result, std::hash<Allele::NucleotideSequence>()(ref_sequence(variant)));
        hash_combine(result, std::hash<Allele::NucleotideSequence>()(alt_sequence(variant)));
        return result;
    }
};

} // namespace octopus

namespace std {
    template <> struct hash<octopus::Variant>
    {
        size_t operator()(const octopus::Variant& variant) const
        {
            return octopus::VariantHash()(variant);
        }
    };
} // namespace std

#endif
