//
//  cigar_scanner.hpp
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__cigar_scanner__
#define __Octopus__cigar_scanner__

#include <vector>
#include <deque>
#include <cstddef>
#include <utility>
#include <functional>

#include "candidate_variant_generator.hpp"
#include "aligned_read.hpp"
#include "variant.hpp"

class ReferenceGenome;
class GenomicRegion;

namespace octopus { namespace core { namespace generators
{
class CigarScanner : public CandidateVariantGenerator
{
public:
    struct Options
    {
        AlignedRead::BaseQuality min_base_quality;
        unsigned min_support = 1;
        Variant::RegionType::Size max_variant_size = 100;
        bool always_include_overlapping_indels = true;
        unsigned max_poor_quality_insertion_bases = 1;
    };
    
    CigarScanner() = delete;
    
    CigarScanner(const ReferenceGenome& reference, Options options);
    
    CigarScanner(const CigarScanner&)            = default;
    CigarScanner& operator=(const CigarScanner&) = default;
    CigarScanner(CigarScanner&&)                 = default;
    CigarScanner& operator=(CigarScanner&&)      = default;
    
    ~CigarScanner() override = default;
    
    bool requires_reads() const noexcept override;
    
    void add_read(const AlignedRead& read) override;
    void add_reads(std::vector<AlignedRead>::const_iterator first,
                   std::vector<AlignedRead>::const_iterator last) override;
    void add_reads(MappableFlatMultiSet<AlignedRead>::const_iterator first,
                   MappableFlatMultiSet<AlignedRead>::const_iterator last) override;
    
    std::vector<Variant> generate_candidates(const GenomicRegion& region) override;
    
    void clear() override;
    
private:
    using NucleotideSequence = AlignedRead::NucleotideSequence;
    using SequenceIterator   = NucleotideSequence::const_iterator;
    using QualitiesIterator  = AlignedRead::BaseQualityVector::const_iterator;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    
    Options options_;
    
    std::function<bool(const Variant&, const Variant&)> match_;
    
    std::deque<Variant> candidates_;
    Variant::RegionType::Size max_seen_candidate_size_;
    
    template <typename T1, typename T2, typename T3>
    void add_candidate(T1&& region, T2&& sequence_removed, T3&& sequence_added);
    void add_snvs_in_match_range(const GenomicRegion& region, SequenceIterator first_base,
                                 SequenceIterator last_base, QualitiesIterator first_quality);
};

template <typename T1, typename T2, typename T3>
void CigarScanner::add_candidate(T1&& region, T2&& sequence_removed,
                                                       T3&& sequence_added)
{
    const auto candidate_size = region_size(region);
    
    if (candidate_size <= options_.max_variant_size) {
        candidates_.emplace_back(std::forward<T1>(region), std::forward<T2>(sequence_removed),
                                 std::forward<T3>(sequence_added));
        
        max_seen_candidate_size_ = std::max(max_seen_candidate_size_, candidate_size);
    }
}
} // namespace generators
} // namespace core
} // namespace octopus

#endif /* defined(__Octopus__cigar_scanner__) */
