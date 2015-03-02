//
//  variant_assembler.cpp
//  Octopus
//
//  Created by Daniel Cooke on 24/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant_assembler.h"

#include "aligned_read.h"
#include "genomic_region.h"
#include "variant.h"

VariantAssembler::VariantAssembler(unsigned k)
:the_assembler_ {k}
{}

void VariantAssembler::add_read(const AlignedRead& a_read)
{
    the_assembler_.add_sequence(a_read.get_sequence(), a_read.get_begin(), Colour::Read);
}

void VariantAssembler::add_reference_sequence(const GenomicRegion& the_region,
                                              const std::string& the_sequence)
{
    the_assembler_.add_sequence(the_sequence, the_region.get_begin(), Colour::Reference);
}

std::set<Variant> VariantAssembler::get_variants(const GenomicRegion& a_region)
{
    std::set<Variant> result {};
    return result;
}

void VariantAssembler::clear() noexcept
{
    the_assembler_.clear();
}