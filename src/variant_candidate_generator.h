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
#include <algorithm> // std::inplace_merge
#include <iterator>  // std::make_move_iterator
#include <cstddef>   // std::size_t
#include <numeric>   // std::accumulate

#include "i_variant_candidate_generator.h"
#include "reference_genome.h"
#include "genomic_region.h"
#include "aligned_read.h"
#include "variant.h"

class VariantCandidateGenerator : public IVariantCandidateGenerator
{
public:
    VariantCandidateGenerator()           = default;
    ~VariantCandidateGenerator() override = default;
    
    VariantCandidateGenerator(const VariantCandidateGenerator&)            = default;
    VariantCandidateGenerator& operator=(const VariantCandidateGenerator&) = default;
    VariantCandidateGenerator(VariantCandidateGenerator&&)                 = default;
    VariantCandidateGenerator& operator=(VariantCandidateGenerator&&)      = default;
    
    void register_generator(std::unique_ptr<IVariantCandidateGenerator> generator);
    void add_read(const AlignedRead& a_read) override;
    void add_reads(ReadIterator first, ReadIterator last) override;
    std::vector<Variant> get_candidates(const GenomicRegion& a_region) override;
    void reserve(std::size_t n) override;
    void clear() override;
    
    std::vector<std::pair<Variant, double>> get_candidates_and_priors(const GenomicRegion& a_region);
    double get_variant_detection_probability(const Variant& a_variant) override;
    
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

inline void VariantCandidateGenerator::add_reads(ReadIterator first, ReadIterator last)
{
    for (auto& generator : generator_list_) {
        generator->add_reads(first, last);
    }
}

inline std::vector<Variant> VariantCandidateGenerator::get_candidates(const GenomicRegion& a_region)
{
    std::vector<Variant> result {};
    
    for (auto& generator : generator_list_) {
        auto generator_result = generator->get_candidates(a_region);
        auto it = result.insert(std::end(result),
                                std::make_move_iterator(std::begin(generator_result)),
                                std::make_move_iterator(std::end(generator_result)));
        std::inplace_merge(std::begin(result), it, std::end(result));
    }
    
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

inline std::vector<std::pair<Variant, double>>
VariantCandidateGenerator::get_candidates_and_priors(const GenomicRegion& a_region)
{
    std::vector<std::pair<Variant, double>> result {};
    
    return result;
}

inline double VariantCandidateGenerator::get_variant_detection_probability(const Variant& a_variant)
{
    double result {1};
    for (auto& generator : generator_list_) {
        result *= generator->get_variant_detection_probability(a_variant);
    }
    return result;
}

#endif /* defined(__Octopus__variant_candidate_generator__) */
