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

#include "i_variant_candidate_generator.h"
#include "variant_factory.h"
#include "aligned_read.h"

class ReferenceGenome;
class GenomicRegion;
class Variant;

class AlignmentCandidateVariantGenerator : public IVariantCandidateGenerator
{
public:
    AlignmentCandidateVariantGenerator() = delete;
    explicit AlignmentCandidateVariantGenerator(ReferenceGenome& the_reference,
                                                VariantFactory& variant_factory,
                                                double generator_confidence);
    ~AlignmentCandidateVariantGenerator() override = default;
    
    AlignmentCandidateVariantGenerator(const AlignmentCandidateVariantGenerator&)            = default;
    AlignmentCandidateVariantGenerator& operator=(const AlignmentCandidateVariantGenerator&) = default;
    AlignmentCandidateVariantGenerator(AlignmentCandidateVariantGenerator&&)                 = default;
    AlignmentCandidateVariantGenerator& operator=(AlignmentCandidateVariantGenerator&&)      = default;
    
    void add_read(const AlignedRead& a_read) override;
    void add_reads(ReadIterator first, ReadIterator last) override;
    std::vector<Variant> get_candidates(const GenomicRegion& a_region) override;
    void reserve(std::size_t n) override;
    void clear() override;
    
    double get_variant_detection_probability(const Variant& a_variant) override;
    
private:
    ReferenceGenome& the_reference_;
    std::vector<Variant> candidates_;
    VariantFactory& variant_factory_;
    double generator_confidence_;
    
    bool is_good_sequence(const std::string& sequence) const noexcept;
    template <typename T1, typename T2, typename T3>
    void add_variant(T1&& the_region, T2&& sequence_removed, T3&& sequence_added);
    void get_variants_in_match_range(const GenomicRegion& the_region,
                                     std::string::const_iterator read_begin,
                                     std::string::const_iterator read_end);
    std::size_t estimate_num_variants(std::size_t num_reads) const noexcept;
};

template <typename T1, typename T2, typename T3>
void AlignmentCandidateVariantGenerator::add_variant(T1&& the_region, T2&& sequence_removed,
                                                     T3&& sequence_added)
{
    candidates_.emplace_back(variant_factory_.make(std::forward<T1>(the_region),
                                                   std::forward<T2>(sequence_removed),
                                                   std::forward<T3>(sequence_added)));
}

#endif /* defined(__Octopus__alignment_candidate_variant_generator__) */
