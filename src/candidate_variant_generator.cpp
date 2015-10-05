//
//  candidate_variant_generator.cpp
//  Octopus
//
//  Created by Daniel Cooke on 14/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "candidate_variant_generator.hpp"

#include <algorithm> // std::inplace_merge

#include "genomic_region.hpp"
#include "aligned_read.hpp"
#include "variant.hpp"
#include "variant_utils.hpp"

namespace Octopus {
    
void CandidateVariantGenerator::register_generator(std::unique_ptr<ICandidateVariantGenerator> generator)
{
    generator_list_.emplace_back(std::move(generator));
}

void CandidateVariantGenerator::add_read(const AlignedRead& a_read)
{
    for (auto& generator : generator_list_) {
        generator->add_read(a_read);
    }
}

void CandidateVariantGenerator::add_reads(std::vector<AlignedRead>::const_iterator first, std::vector<AlignedRead>::const_iterator last)
{
    for (auto& generator : generator_list_) {
        generator->add_reads(first, last);
    }
}

void CandidateVariantGenerator::add_reads(MappableSet<AlignedRead>::const_iterator first, MappableSet<AlignedRead>::const_iterator last)
{
    for (auto& generator : generator_list_) {
        generator->add_reads(first, last);
    }
}

std::vector<Variant> CandidateVariantGenerator::get_candidates(const GenomicRegion& a_region)
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

void CandidateVariantGenerator::reserve(size_t n)
{
    for (auto& generator : generator_list_) {
        generator->reserve(n);
    }
}

void CandidateVariantGenerator::clear()
{
    for (auto& generator : generator_list_) {
        generator->clear();
    }
}

} // namespace Octopus