//
//  variant_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 15/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant_caller.h"

#include "genomic_region.h"
#include "mappable.h"
#include "read_utils.h"
#include "vcf_record.h"

#include <iostream> // TEST

VariantCaller::VariantCaller(ReferenceGenome& reference, ReadManager& read_manager,
                             ReadFilter read_filter, ReadTransform read_transform,
                             CandidateVariantGenerator& candidate_generator)
:
reference_ {reference},
read_manager_ {read_manager},
read_filter_ {read_filter},
read_transform_ {read_transform},
candidate_generator_ {candidate_generator}
{}

bool VariantCaller::done_calling(const GenomicRegion& region) const noexcept
{
    return empty(region);
}

std::vector<VcfRecord> VariantCaller::call_variants(const GenomicRegion& region)
{
    auto current_region = get_init_region(region);
    
    std::vector<VcfRecord> result {};
    
    while (!done_calling(current_region)) {
        auto reads = make_mappable_map(read_manager_.fetch_reads(current_region));
        
        auto good_reads = filter_reads(reads, read_filter_).first;
        
        //std::cout << "found " << good_reads[samples.front()].size() << " good reads" << std::endl;
        
        transform_reads(good_reads, read_transform_);
        
        std::cout << "transformed reads" << std::endl;
        
        add_reads(good_reads, candidate_generator_);
        
        auto candidates = candidate_generator_.get_candidates(current_region);
        
        candidate_generator_.clear();
        
        std::cout << "found " << candidates.size() << " candidates" << std::endl;
        
        auto new_variants = call_variants(current_region, candidates, good_reads);
        
        result.insert(std::end(result), std::begin(new_variants), std::end(new_variants));
        
        current_region = get_next_region(current_region);
    }
    
    return result;
}
