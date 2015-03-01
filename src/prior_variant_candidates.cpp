//
//  PriorVariantCandidates.cpp
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "prior_variant_candidates.h"
#include "variant.h"
#include "variant_file.h"

PriorVariantCandidates::PriorVariantCandidates(VariantFile& a_variant_source)
: a_variant_file_ {a_variant_source}
{}

std::set<Variant> PriorVariantCandidates::get_candidates(const GenomicRegion& a_region)
{
    return a_variant_file_.fetch_variants(a_region);
}
