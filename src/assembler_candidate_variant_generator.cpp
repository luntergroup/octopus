//
//  assembler_candidate_variant_generator.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "assembler_candidate_variant_generator.hpp"

#include <algorithm>

#include "reference_genome.hpp"
#include "aligned_read.hpp"
#include "variant.hpp"

namespace Octopus {
    
AssemblerCandidateVariantGenerator::AssemblerCandidateVariantGenerator(ReferenceGenome& reference,
                                                                       unsigned kmer_size,
                                                                       SizeType max_variant_size)
:
reference_ {reference},
assembler_ {kmer_size},
max_variant_size_ {max_variant_size}
{}

void AssemblerCandidateVariantGenerator::add_read(const AlignedRead& read)
{
    assembler_.add_read(read);
}

void AssemblerCandidateVariantGenerator::add_reads(std::vector<AlignedRead>::const_iterator first,
                                                   std::vector<AlignedRead>::const_iterator last)
{
    std::for_each(first, last, [this] (const auto& read ) { add_read(read); });
}

void AssemblerCandidateVariantGenerator::add_reads(MappableSet<AlignedRead>::const_iterator first,
                                                   MappableSet<AlignedRead>::const_iterator last)
{
    std::for_each(first, last, [this] (const auto& read ) { add_read(read); });
}

std::vector<Variant> AssemblerCandidateVariantGenerator::get_candidates(const GenomicRegion& region)
{
    auto reference_sequence = reference_.get_sequence(region);
    assembler_.add_reference_sequence(region, reference_sequence);
    return assembler_.get_variants(region);
}

void AssemblerCandidateVariantGenerator::clear()
{
    assembler_.clear();
}

} // namespace Octopus
