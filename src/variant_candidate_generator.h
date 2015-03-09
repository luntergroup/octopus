//
//  variant_candidate_generator.h
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__variant_candidate_generator__
#define __Octopus__variant_candidate_generator__

#include <vector>
#include <memory>    // std::unique_ptr
#include <algorithm> // std::for_each, std::sort
#include <iterator>  // std::make_move_iterator
#include <cstddef>   // std::size_t

#include "i_variant_candidate_generator.h"
#include "reference_genome.h"
#include "genomic_region.h"
#include "aligned_read.h"
#include "variant.h"

class VariantCandidateGenerator : public IVariantCandidateGenerator
{
public:
    VariantCandidateGenerator() = default;
    ~VariantCandidateGenerator() override = default;
    
    VariantCandidateGenerator(const VariantCandidateGenerator&)            = default;
    VariantCandidateGenerator& operator=(const VariantCandidateGenerator&) = default;
    VariantCandidateGenerator(VariantCandidateGenerator&&)                 = default;
    VariantCandidateGenerator& operator=(VariantCandidateGenerator&&)      = default;
    
    void register_generator(std::unique_ptr<IVariantCandidateGenerator> generator);
    void add_read(const AlignedRead& a_read) override;
    std::vector<Variant> get_candidates(const GenomicRegion& a_region) override;
    void reserve(std::size_t n) override;
    void clear() override;
    template <typename ReadIterator> void add_reads(ReadIterator first, ReadIterator last);
    
private:
    std::vector<std::unique_ptr<IVariantCandidateGenerator>> generator_list_;
};

inline void VariantCandidateGenerator::register_generator(std::unique_ptr<IVariantCandidateGenerator> generator)
{
    generator_list_.emplace_back(std::move(generator));
}

inline void VariantCandidateGenerator::add_read(const AlignedRead& a_read)
{
    for (auto& generator : generator_list_) {
        generator->add_read(a_read);
    }
}

inline std::vector<Variant> VariantCandidateGenerator::get_candidates(const GenomicRegion& a_region)
{
    std::vector<Variant> result {};
    for (auto& generator : generator_list_) {
        auto generator_result = generator->get_candidates(a_region);
        result.insert(std::end(result),
                      std::make_move_iterator(std::begin(generator_result)),
                      std::make_move_iterator(std::end(generator_result))
                      );
    }
    std::sort(std::begin(result), std::end(result));
    return result;
}

inline void VariantCandidateGenerator::reserve(std::size_t n)
{
    for (auto& generator : generator_list_) {
        generator->reserve(n);
    }
}

inline void VariantCandidateGenerator::clear()
{
    for (auto& generator : generator_list_) {
        generator->clear();
    }
}

template <typename ReadIterator>
void VariantCandidateGenerator::add_reads(ReadIterator first, ReadIterator last)
{
    std::for_each(first, last, [this] (const auto& a_read ) { add_read(a_read); });
}

#endif /* defined(__Octopus__variant_candidate_generator__) */
