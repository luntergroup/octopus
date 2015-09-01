//
//  assembler_candidate_variant_generator.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "assembler_candidate_variant_generator.h"

#include <algorithm> // std::for_each

#include "reference_genome.h"
#include "aligned_read.h"
#include "variant.h"

AssemblerCandidateVariantGenerator::AssemblerCandidateVariantGenerator(ReferenceGenome& reference,
                                                                       unsigned kmer_size,
                                                                       double generator_confidence)
:
reference_ {reference},
the_variant_assembler_ {kmer_size},
generator_confidence_ {generator_confidence}
{}

void AssemblerCandidateVariantGenerator::add_read(const AlignedRead& a_read)
{
    the_variant_assembler_.add_read(a_read);
}

void AssemblerCandidateVariantGenerator::add_reads(std::vector<AlignedRead>::const_iterator first, std::vector<AlignedRead>::const_iterator last)
{
    std::for_each(first, last, [this] (const auto& a_read ) { add_read(a_read); });
}

void AssemblerCandidateVariantGenerator::add_reads(MappableSet<AlignedRead>::const_iterator first, MappableSet<AlignedRead>::const_iterator last)
{
    std::for_each(first, last, [this] (const auto& a_read ) { add_read(a_read); });
}

std::vector<Variant> AssemblerCandidateVariantGenerator::get_candidates(const GenomicRegion& region)
{
    auto reference_sequence = reference_.get_sequence(region);
    the_variant_assembler_.add_reference_sequence(region, reference_sequence);
    return the_variant_assembler_.get_variants(region);
}

void AssemblerCandidateVariantGenerator::clear()
{
    the_variant_assembler_.clear();
}
