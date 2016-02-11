//
//  variant_assembler.cpp
//  Octopus
//
//  Created by Daniel Cooke on 24/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant_assembler.hpp"

#include "aligned_read.hpp"
#include "genomic_region.hpp"
#include "variant.hpp"

VariantAssembler::VariantAssembler(unsigned k)
:
de_bruijn_graph_ {k}
{}

void VariantAssembler::add_read(const AlignedRead& read)
{
    de_bruijn_graph_.add_sequence(read.get_sequence(), region_begin(read), Colour::Read);
}

void VariantAssembler::add_reference_sequence(const GenomicRegion& region, const std::string& sequence)
{
    de_bruijn_graph_.add_sequence(sequence, region_begin(region), Colour::Reference);
}

std::vector<Variant> VariantAssembler::get_variants(const GenomicRegion& region)
{
    std::vector<Variant> result {};
    
    de_bruijn_graph_.print_kmers(9378586, Colour::Reference);
    
    de_bruijn_graph_.get_contigs(1);
    
    return result;
}

void VariantAssembler::clear() noexcept
{
    de_bruijn_graph_.clear();
}