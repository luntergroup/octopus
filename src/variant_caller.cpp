//
//  variant_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 15/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant_caller.hpp"

#include "genomic_region.hpp"
#include "mappable.hpp"
#include "read_utils.hpp"
#include "variant_utils.hpp"
#include "vcf_record.hpp"

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
        
        transform_reads(good_reads, read_transform_);
        
        add_reads(good_reads, candidate_generator_);
        
        auto candidates = unique_left_align(candidate_generator_.get_candidates(current_region), reference_);
        
        candidate_generator_.clear();
        
        //std::cout << "found " << candidates.size() << " candidates" << std::endl;
        
        auto calls_in_region = call_variants(current_region, candidates, good_reads);
        
        result.insert(std::end(result), std::begin(calls_in_region), std::end(calls_in_region));
        
        current_region = get_next_region(current_region);
    }
    
    return result;
}
