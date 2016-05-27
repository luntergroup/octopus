//
//  alignment_candidate_variant_generator.hpp
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__alignment_candidate_variant_generator__
#define __Octopus__alignment_candidate_variant_generator__

#include <vector>
#include <deque>
#include <cstddef>
#include <utility>
#include <functional>

#include "i_candidate_variant_generator.hpp"
#include "aligned_read.hpp"

class ReferenceGenome;
class GenomicRegion;
class Variant;

namespace Octopus {
    
class AlignmentCandidateVariantGenerator : public ICandidateVariantGenerator
{
public:
    using QualityType = AlignedRead::QualityType;
    using SizeType    = GenomicRegion::SizeType;
    
    AlignmentCandidateVariantGenerator() = delete;
    
    AlignmentCandidateVariantGenerator(const ReferenceGenome& reference,
                                       QualityType min_base_quality = 0,
                                       unsigned min_support = 2,
                                       SizeType max_variant_size = 100);
    
    ~AlignmentCandidateVariantGenerator() override = default;
    
    AlignmentCandidateVariantGenerator(const AlignmentCandidateVariantGenerator&)            = default;
    AlignmentCandidateVariantGenerator& operator=(const AlignmentCandidateVariantGenerator&) = default;
    AlignmentCandidateVariantGenerator(AlignmentCandidateVariantGenerator&&)                 = default;
    AlignmentCandidateVariantGenerator& operator=(AlignmentCandidateVariantGenerator&&)      = default;
    
    bool requires_reads() const noexcept override;
    void add_read(const AlignedRead& read) override;
    void add_reads(std::vector<AlignedRead>::const_iterator first,
                   std::vector<AlignedRead>::const_iterator last) override;
    void add_reads(MappableFlatMultiSet<AlignedRead>::const_iterator first,
                   MappableFlatMultiSet<AlignedRead>::const_iterator last) override;
    std::vector<Variant> generate_candidates(const GenomicRegion& region) override;
    void clear() override;
    
private:
    using SequenceType      = AlignedRead::SequenceType;
    using SequenceIterator  = SequenceType::const_iterator;
    using QualitiesIterator = AlignedRead::Qualities::const_iterator;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    QualityType min_base_quality_;
    unsigned max_poor_quality_insertion_bases_;
    unsigned min_support_;
    SizeType max_variant_size_;
    
    std::deque<Variant> candidates_;
    SizeType max_seen_candidate_size_;
    
    template <typename T1, typename T2, typename T3>
    void add_candidate(T1&& region, T2&& sequence_removed, T3&& sequence_added);
    void add_snvs_in_match_range(const GenomicRegion& region, SequenceIterator first_base,
                                 SequenceIterator last_base, QualitiesIterator first_quality);
};

template <typename T1, typename T2, typename T3>
void AlignmentCandidateVariantGenerator::add_candidate(T1&& region, T2&& sequence_removed,
                                                       T3&& sequence_added)
{
    const auto candidate_size = region_size(region);
    
    if (candidate_size <= max_variant_size_) {
        candidates_.emplace_back(std::forward<T1>(region), std::forward<T2>(sequence_removed),
                                 std::forward<T3>(sequence_added));
        
        max_seen_candidate_size_ = std::max(max_seen_candidate_size_, candidate_size);
    }
}

} // namespace Octopus

#endif /* defined(__Octopus__alignment_candidate_variant_generator__) */
