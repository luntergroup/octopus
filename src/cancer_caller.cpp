//
//  cancer_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 16/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "cancer_caller.h"

#include <unordered_map>
#include <numeric>

#include "genomic_region.h"
#include "read_manager.h"
#include "allele.h"
#include "variant.h"
#include "haplotype.h"
#include "genotype.h"
#include "haplotype_tree.h"
#include "search_regions.h"
#include "vcf_record.h"

#include "genotype_model.h"
#include "population_genotype_model.h"

#include <iostream> // TEST

CancerVariantCaller::CancerVariantCaller(ReferenceGenome& reference, ReadManager& read_manager,
                                         ReadFilter read_filter, ReadTransform read_transform,
                                         CandidateVariantGenerator& candidate_generator)
:
VariantCaller {reference, read_manager, read_filter, read_transform, candidate_generator}
{}

GenomicRegion CancerVariantCaller::get_init_region(const GenomicRegion& region)
{
    return region;
}

GenomicRegion CancerVariantCaller::get_next_region(const GenomicRegion& current_region)
{
    return GenomicRegion {"TEST", 0, 0};
}

std::vector<VcfRecord> CancerVariantCaller::call_variants(const GenomicRegion& region,
                                                          const std::vector<Variant>& candidates,
                                                          const ReadMap& reads)
{
    std::vector<VcfRecord> result {};
    
    return result;
}
