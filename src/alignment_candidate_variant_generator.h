//
//  alignment_candidate_variant_generator.h
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__alignment_candidate_variant_generator__
#define __Octopus__alignment_candidate_variant_generator__

#include <vector>
#include <cstddef> // std::size_t
#include <algorithm> // std::for_each, std::sort, std::unique, std::lower_bound, std::upper_bound, std::max

#include "i_candidate_variant_generator.h"
#include "aligned_read.h"

class ReferenceGenome;
class GenomicRegion;
class Variant;

class AlignmentCandidateVariantGenerator : public ICandidateVariantGenerator
{
public:
    using QualityType = AlignedRead::QualityType;
    using SizeType    = GenomicRegion::SizeType;
    
    AlignmentCandidateVariantGenerator() = delete;
    explicit AlignmentCandidateVariantGenerator(ReferenceGenome& reference,
                                                QualityType min_base_quality = 0,
                                                SizeType max_variant_size = 100);
    ~AlignmentCandidateVariantGenerator() override = default;
    
    AlignmentCandidateVariantGenerator(const AlignmentCandidateVariantGenerator&)            = default;
    AlignmentCandidateVariantGenerator& operator=(const AlignmentCandidateVariantGenerator&) = default;
    AlignmentCandidateVariantGenerator(AlignmentCandidateVariantGenerator&&)                 = default;
    AlignmentCandidateVariantGenerator& operator=(AlignmentCandidateVariantGenerator&&)      = default;
    
    void add_read(const AlignedRead& read) override;
    void add_reads(std::vector<AlignedRead>::const_iterator first, std::vector<AlignedRead>::const_iterator last) override;
    void add_reads(MappableSet<AlignedRead>::const_iterator first, MappableSet<AlignedRead>::const_iterator last) override;
    std::vector<Variant> get_candidates(const GenomicRegion& region) override;
    void reserve(std::size_t n) override;
    void clear() override;
    
private:
    using SequenceType      = AlignedRead::SequenceType;
    using SequenceIterator  = SequenceType::const_iterator;
    using QualitiesIterator = AlignedRead::Qualities::const_iterator;
    
    ReferenceGenome& reference_;
    QualityType min_base_quality_;
    SizeType max_variant_size_;
    
    std::vector<Variant> candidates_;
    bool are_candidates_sorted_;
    SizeType max_seen_candidate_size_;
    
    bool is_good_sequence(const SequenceType& sequence) const noexcept;
    template <typename T1, typename T2, typename T3>
    void add_variant(T1&& the_region, T2&& sequence_removed, T3&& sequence_added);
    void get_snvs_in_match_range(const GenomicRegion& the_region, SequenceIterator first_base,
                                 SequenceIterator last_base, QualitiesIterator first_quality);
    std::size_t estimate_num_variants(std::size_t num_reads) const noexcept;
};

template <typename T1, typename T2, typename T3>
void AlignmentCandidateVariantGenerator::add_variant(T1&& the_region, T2&& sequence_removed,
                                                     T3&& sequence_added)
{
    auto candidate_size = size(the_region);
    if (candidate_size <= max_variant_size_) {
        candidates_.emplace_back(std::forward<T1>(the_region), std::forward<T2>(sequence_removed),
                                 std::forward<T3>(sequence_added));
        max_seen_candidate_size_ = std::max(max_seen_candidate_size_, candidate_size);
        are_candidates_sorted_ = false;
    }
}

#endif /* defined(__Octopus__alignment_candidate_variant_generator__) */
