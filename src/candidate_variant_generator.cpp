//
//  candidate_variant_generator.cpp
//  Octopus
//
//  Created by Daniel Cooke on 14/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "candidate_variant_generator.hpp"

#include <algorithm>
#include <iterator>

#include "genomic_region.hpp"
#include "aligned_read.hpp"
#include "variant.hpp"

namespace Octopus {
    
void CandidateVariantGenerator::register_generator(std::unique_ptr<ICandidateVariantGenerator> generator)
{
    generators_.emplace_back(std::move(generator));
}

bool CandidateVariantGenerator::requires_reads() const noexcept
{
    return std::any_of(std::cbegin(generators_), std::cend(generators_),
                       [] (const auto& generator) { return generator->requires_reads(); });
}

void CandidateVariantGenerator::add_read(const AlignedRead& read)
{
    for (auto& generator : generators_) {
        generator->add_read(read);
    }
}

void CandidateVariantGenerator::add_reads(std::vector<AlignedRead>::const_iterator first,
                                          std::vector<AlignedRead>::const_iterator last)
{
    for (auto& generator : generators_) {
        generator->add_reads(first, last);
    }
}

void CandidateVariantGenerator::add_reads(MappableFlatMultiSet<AlignedRead>::const_iterator first,
                                          MappableFlatMultiSet<AlignedRead>::const_iterator last)
{
    for (auto& generator : generators_) {
        generator->add_reads(first, last);
    }
}

std::vector<Variant> CandidateVariantGenerator::generate_candidates(const GenomicRegion& region)
{
    std::vector<Variant> result {};
    
    for (auto& generator : generators_) {
        auto generator_result = generator->generate_candidates(region); // results are sorted
        const auto it = result.insert(std::end(result),
                                      std::make_move_iterator(std::begin(generator_result)),
                                      std::make_move_iterator(std::end(generator_result)));
        std::inplace_merge(std::begin(result), it, std::end(result));
    }
    
    remove_duplicates(result);
    
    return result;
}

void CandidateVariantGenerator::reserve(const size_t n)
{
    for (auto& generator : generators_) {
        generator->reserve(n);
    }
}

void CandidateVariantGenerator::clear()
{
    for (auto& generator : generators_) {
        generator->clear();
    }
}

} // namespace Octopus