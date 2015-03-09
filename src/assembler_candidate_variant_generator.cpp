//
//  assembler_candidate_variant_generator.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "assembler_candidate_variant_generator.h"
#include "reference_genome.h"
#include "aligned_read.h"
#include "variant.h"

AssemblerCandidateVariantGenerator::AssemblerCandidateVariantGenerator(ReferenceGenome& the_reference,
                                                                       unsigned kmer_size)
:the_reference_ {the_reference},
the_variant_assembler_ {kmer_size}
{}

void AssemblerCandidateVariantGenerator::add_read(const AlignedRead& a_read)
{
    the_variant_assembler_.add_read(a_read);
}

std::vector<Variant> AssemblerCandidateVariantGenerator::get_candidates(const GenomicRegion& a_region)
{
    auto reference_sequence = the_reference_.get_sequence(a_region);
    the_variant_assembler_.add_reference_sequence(a_region, reference_sequence);
    return the_variant_assembler_.get_variants(a_region);
}

void AssemblerCandidateVariantGenerator::reserve(std::size_t n) {}

void AssemblerCandidateVariantGenerator::clear()
{
    the_variant_assembler_.clear();
}
