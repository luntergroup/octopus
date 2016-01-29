//
//  allele.hpp
//  Octopus
//
//  Created by Daniel Cooke on 10/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_allele_hpp
#define Octopus_allele_hpp

#include <string>
#include <ostream>
#include <utility>

#include <boost/functional/hash.hpp>

#include "genomic_region.hpp"
#include "comparable.hpp"
#include "mappable.hpp"
#include "reference_genome.hpp"

class Allele : public Comparable<Allele>, public Mappable<Allele>
{
public:
    using SizeType     = GenomicRegion::SizeType;
    using SequenceType = ReferenceGenome::SequenceType;
    
    Allele() = default;
    template <typename GenomicRegion_, typename SequenceType_>
    Allele(GenomicRegion_&& reference_region, SequenceType_&& sequence);
    template <typename StringType_, typename SequenceType_>
    Allele(StringType_&& contig_name, SizeType begin_pos, SequenceType_&& sequence);
    ~Allele() = default;
    
    Allele(const Allele&)            = default;
    Allele& operator=(const Allele&) = default;
    Allele(Allele&&)                 = default;
    Allele& operator=(Allele&&)      = default;
    
    const GenomicRegion& get_region() const noexcept;
    const SequenceType& get_sequence() const noexcept;
    
private:
    SequenceType sequence_;
    GenomicRegion reference_region_;
};

template <typename GenomicRegion_, typename SequenceType_>
Allele::Allele(GenomicRegion_&& reference_region, SequenceType_&& sequence)
:
sequence_ {std::forward<SequenceType_>(sequence)},
reference_region_ {std::forward<GenomicRegion_>(reference_region)}
{}

template <typename StringType_, typename SequenceType_>
Allele::Allele(StringType_&& contig_name, SizeType begin_pos, SequenceType_&& sequence)
:
sequence_ {std::forward<SequenceType_>(sequence)},
reference_region_ {std::forward<StringType_>(contig_name), begin_pos,
        static_cast<SizeType>(begin_pos + sequence_.size())}
{}

// non-member methods

Allele::SizeType sequence_size(const Allele& allele) noexcept;

bool is_reference(const Allele& allele, const ReferenceGenome& reference);

Allele get_reference_allele(const GenomicRegion& region, const ReferenceGenome& reference);

std::vector<Allele> get_reference_alleles(const std::vector<GenomicRegion>& regions,
                                          const ReferenceGenome& reference);

std::vector<Allele> get_positional_reference_alleles(const GenomicRegion& region,
                                                     const ReferenceGenome& reference);

bool contains(const Allele& lhs, const Allele& rhs);

Allele splice(const Allele& allele, const GenomicRegion& region);

bool is_insertion(const Allele& allele);
bool is_deletion(const Allele& allele);
bool is_indel(const Allele& allele);

std::vector<Allele> decompose(const Allele& allele);

bool operator==(const Allele& lhs, const Allele& rhs);
bool operator<(const Allele& lhs, const Allele& rhs);

namespace std {
    template <> struct hash<Allele>
    {
        size_t operator()(const Allele& allele) const
        {
            using boost::hash_combine;
            size_t result {};
            hash_combine(result, hash<GenomicRegion>()(allele.get_region()));
            hash_combine(result, hash<Allele::SequenceType>()(allele.get_sequence()));
            return result;
        }
    };
} // namespace std

namespace boost
{
    template <> struct hash<Allele> : std::unary_function<Allele, std::size_t>
    {
        std::size_t operator()(const Allele& a) const
        {
            return std::hash<Allele>()(a);
        }
    };
} // namespace boost

std::ostream& operator<<(std::ostream& os, const Allele& allele);

#endif
