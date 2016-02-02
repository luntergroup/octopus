//
//  contig_allele.hpp
//  Octopus
//
//  Created by Daniel Cooke on 21/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef contig_allele_hpp
#define contig_allele_hpp

#include <string>
#include <cstddef>
#include <ostream>

#include "comparable.hpp"
#include "contig_region.hpp"
#include "genomic_region.hpp"
#include "mappable.hpp"
#include "allele.hpp"

/*
 ContigAllele is like Allele but uses ContigRegion rather than GenomicRegion, this makes it much
 less flexible, but more efficient to use within other classes as a set with a common contig.
 */
class ContigAllele : public Comparable<ContigAllele>, public Mappable<ContigAllele>
{
public:
    using SizeType     = ContigRegion::SizeType;
    using SequenceType = Allele::SequenceType;
    
    ContigAllele() = default;
    template <typename SequenceType_> explicit ContigAllele(ContigRegion region, SequenceType_&& sequence);
    template <typename SequenceType_> explicit ContigAllele(const GenomicRegion& region, SequenceType_&& sequence);
    explicit ContigAllele(const Allele& allele);
    ~ContigAllele() = default;
    
    ContigAllele(const ContigAllele&)            = default;
    ContigAllele& operator=(const ContigAllele&) = default;
    ContigAllele(ContigAllele&&)                 = default;
    ContigAllele& operator=(ContigAllele&&)      = default;
    
    const ContigRegion& get_region() const noexcept;
    const SequenceType& get_sequence() const noexcept;
    
private:
    SequenceType sequence_;
    ContigRegion region_;
};

template <typename SequenceType_>
ContigAllele::ContigAllele(ContigRegion region, SequenceType_&& sequence)
:
region_ {region},
sequence_ {std::forward<SequenceType_>(sequence)}
{}

template <typename SequenceType_>
ContigAllele::ContigAllele(const GenomicRegion& region, SequenceType_&& sequence)
:
region_ {region.get_contig_region()},
sequence_ {std::forward<SequenceType_>(sequence)}
{}

// non-member methods

ContigAllele::SizeType sequence_size(const ContigAllele& allele) noexcept;

bool contains(const ContigAllele& lhs, const ContigAllele& rhs);

ContigAllele splice(const ContigAllele& allele, const ContigRegion& region);

bool is_insertion(const ContigAllele& allele);
bool is_deletion(const ContigAllele& allele);
bool is_indel(const ContigAllele& allele);

bool operator==(const ContigAllele& lhs, const ContigAllele& rhs);
bool operator<(const ContigAllele& lhs, const ContigAllele& rhs);

namespace std {
    template <> struct hash<ContigAllele>
    {
        size_t operator()(const ContigAllele& allele) const
        {
            using boost::hash_combine;
            size_t result {0};
            hash_combine(result, hash<ContigRegion>()(allele.get_region()));
            hash_combine(result, hash<ContigAllele::SequenceType>()(allele.get_sequence()));
            return result;
        }
    };
} // namespace std

namespace boost
{
    template <> struct hash<ContigAllele> : std::unary_function<ContigAllele, std::size_t>
    {
        std::size_t operator()(const ContigAllele& a) const
        {
            return std::hash<ContigAllele>()(a);
        }
    };
} // namespace boost

std::ostream& operator<<(std::ostream& os, const ContigAllele& allele);

#endif /* contig_allele_hpp */
