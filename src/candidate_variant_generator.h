//
//  candidate_variant_generator.h
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__candidate_variant_generator__
#define __Octopus__candidate_variant_generator__

#include <vector>
#include <memory>    // std::unique_ptr
#include <algorithm> // std::inplace_merge
#include <iterator>  // std::make_move_iterator
#include <cstddef>   // std::size_t
#include <numeric>   // std::accumulate

#include "i_candidate_variant_generator.h"
#include "reference_genome.h"
#include "genomic_region.h"
#include "aligned_read.h"
#include "variant.h"
#include "variant_utils.h"

class CandidateVariantGenerator : public ICandidateVariantGenerator
{
public:
    CandidateVariantGenerator()           = default;
    ~CandidateVariantGenerator() override = default;
    
    CandidateVariantGenerator(const CandidateVariantGenerator&)            = default;
    CandidateVariantGenerator& operator=(const CandidateVariantGenerator&) = default;
    CandidateVariantGenerator(CandidateVariantGenerator&&)                 = default;
    CandidateVariantGenerator& operator=(CandidateVariantGenerator&&)      = default;
    
    void register_generator(std::unique_ptr<ICandidateVariantGenerator> generator);
    void add_read(const AlignedRead& a_read) override;
    void add_reads(std::vector<AlignedRead>::const_iterator first, std::vector<AlignedRead>::const_iterator last) override;
    void add_reads(MappableSet<AlignedRead>::const_iterator first, MappableSet<AlignedRead>::const_iterator last) override;
    std::vector<Variant> get_candidates(const GenomicRegion& a_region) override;
    void reserve(std::size_t n) override;
    void clear() override;
    
private:
    std::vector<std::unique_ptr<ICandidateVariantGenerator>> generator_list_;
};

inline void CandidateVariantGenerator::register_generator(std::unique_ptr<ICandidateVariantGenerator> generator)
{
    generator_list_.emplace_back(std::move(generator));
}

inline void CandidateVariantGenerator::add_read(const AlignedRead& a_read)
{
    for (auto& generator : generator_list_) {
        generator->add_read(a_read);
    }
}

inline void CandidateVariantGenerator::add_reads(std::vector<AlignedRead>::const_iterator first, std::vector<AlignedRead>::const_iterator last)
{
    for (auto& generator : generator_list_) {
        generator->add_reads(first, last);
    }
}

inline void CandidateVariantGenerator::add_reads(MappableSet<AlignedRead>::const_iterator first, MappableSet<AlignedRead>::const_iterator last)
{
    for (auto& generator : generator_list_) {
        generator->add_reads(first, last);
    }
}

inline std::vector<Variant> CandidateVariantGenerator::get_candidates(const GenomicRegion& a_region)
{
    std::vector<Variant> result {};
    
    for (auto& generator : generator_list_) {
        auto generator_result = generator->get_candidates(a_region);
        auto it = result.insert(std::end(result),
                                std::make_move_iterator(std::begin(generator_result)),
                                std::make_move_iterator(std::end(generator_result)));
        std::inplace_merge(std::begin(result), it, std::end(result));
    }
    
    remove_duplicates(result);
    
    return result;
}

inline void CandidateVariantGenerator::reserve(std::size_t n)
{
    for (auto& generator : generator_list_) {
        generator->reserve(n);
    }
}

inline void CandidateVariantGenerator::clear()
{
    for (auto& generator : generator_list_) {
        generator->clear();
    }
}

#endif /* defined(__Octopus__variant_candidate_generator__) */
