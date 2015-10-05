//
//  variant_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 15/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant_caller.hpp"

#include <algorithm> // std::stable_sort

#include "genomic_region.hpp"
#include "mappable.hpp"
#include "read_utils.hpp"
#include "variant_utils.hpp"
#include "vcf_record.hpp"

#include "mappable_algorithms.hpp"

#include <iostream> // TEST

namespace Octopus
{

    // public methods
    
    VariantCaller::VariantCaller(ReferenceGenome& reference, CandidateVariantGenerator& candidate_generator,
                                 RefCallType refcall_type)
    :
    reference_ {reference},
    candidate_generator_ {candidate_generator},
    refcall_type_ {refcall_type}
    {}

    std::string VariantCaller::get_details() const
    {
        return do_get_details();
    }
    
    std::vector<VcfRecord> VariantCaller::call_variants(const GenomicRegion& region, ReadMap reads)
    {
        add_reads(reads, candidate_generator_);
        
        auto candidates = unique_left_align(candidate_generator_.get_candidates(region), reference_);
        
        candidate_generator_.clear();
        
        std::cout << "found " << candidates.size() << " candidates" << std::endl;
        
        auto current_region = get_init_region(region, reads, candidates);
        
        std::vector<VcfRecord> result {};
        
        while (!done_calling(current_region)) {
            std::cout << "processing sub-region " << current_region << std::endl;
            
            auto overlapped = overlap_range(std::cbegin(candidates), std::cend(candidates), current_region);
            
            std::vector<Variant> sub_candidates {std::begin(overlapped), std::end(overlapped)};
            
            auto reads_in_region = copy_overlapped(reads, current_region);
            
            auto calls_in_region = call_variants(current_region, sub_candidates, reads_in_region);
            
            result.insert(std::end(result), std::make_move_iterator(std::begin(calls_in_region)),
                          std::make_move_iterator(std::end(calls_in_region)));
            
            current_region = get_next_region(current_region, reads, candidates);
        }
        
        return result;
    }

    // protected methods

    bool VariantCaller::refcalls_requested() const noexcept
    {
        return refcall_type_ != RefCallType::None;
    }

    // private methods

    bool VariantCaller::done_calling(const GenomicRegion& region) const noexcept
    {
        return empty(region);
    }
    
    // non-member methods
    
    std::vector<Allele>
    generate_callable_alleles(const GenomicRegion& region, const std::vector<Variant>& variants,
                              VariantCaller::RefCallType refcall_type, ReferenceGenome& reference)
    {
        using std::begin; using std::end; using std::make_move_iterator; using std::back_inserter;
        
        if (empty(region)) return {};
        
        if (variants.empty()) {
            switch (refcall_type) {
                case VariantCaller::RefCallType::Positional:
                    return get_positional_reference_alleles(region, reference);
                case VariantCaller::RefCallType::Blocked:
                    return {get_reference_allele(region, reference)};
                default:
                    return {};
            }
        }
        
        auto variant_alleles = decompose(variants);
        
        if (refcall_type == VariantCaller::RefCallType::None) return variant_alleles;
        
        auto covered_regions   = get_covered_regions(variants);
        auto uncovered_regions = get_all_intervening(covered_regions, region);
        
        std::vector<Allele> result {};
        
        if (refcall_type == VariantCaller::RefCallType::Blocked) {
            auto reference_alleles = get_reference_alleles(uncovered_regions, reference);
            result.reserve(reference_alleles.size() + variant_alleles.size());
            std::merge(make_move_iterator(begin(reference_alleles)),make_move_iterator(end(reference_alleles)),
                       make_move_iterator(begin(variant_alleles)), make_move_iterator(end(variant_alleles)),
                       back_inserter(result));
        } else {
            result.reserve(variant_alleles.size() + sum_sizes(uncovered_regions));
            
            auto uncovered_itr = begin(uncovered_regions);
            auto uncovered_end = end(uncovered_regions);
            
            for (auto&& variant_allele : variant_alleles) {
                if (uncovered_itr != uncovered_end && begins_before(*uncovered_itr, variant_allele)) {
                    auto alleles = get_positional_reference_alleles(*uncovered_itr, reference);
                    result.insert(end(result), make_move_iterator(begin(alleles)),
                                  make_move_iterator(end(alleles)));
                    std::advance(uncovered_itr, 1);
                }
                result.push_back(std::move(variant_allele));
            }
            
            if (uncovered_itr != uncovered_end) {
                auto alleles = get_positional_reference_alleles(*uncovered_itr, reference);
                result.insert(end(result), make_move_iterator(begin(alleles)),
                              make_move_iterator(end(alleles)));
            }
        }
        
//        std::stable_sort(begin(result), end(result), [] (const auto& lhs, const auto& rhs) {
//            return get_region(lhs) < get_region(rhs);
//        });
        
        return result;
    }
    
} // namespace Octopus
