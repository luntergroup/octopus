//
//  variant_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 15/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant_caller.hpp"

#include <algorithm> // std::merge, std::copy

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
    
    size_t VariantCaller::num_buffered_reads() const noexcept
    {
        return 0;
    }
    
    std::vector<VcfRecord> VariantCaller::call_variants(const GenomicRegion& region, ReadMap reads)
    {
        add_reads(reads, candidate_generator_);
        
        auto candidates = unique_left_align(candidate_generator_.get_candidates(region), reference_);
        
        candidate_generator_.clear();
        
        std::cout << "found " << candidates.size() << " candidates in " << count_reads(reads) << " reads" << std::endl;
        
        std::cout << "candidates are:" << std::endl;
        for (const auto& c : candidates) std::cout << c << std::endl;
        
        return call_variants(region, candidates, reads);
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
        
        return result;
    }
    
    void append_allele(std::vector<Allele>& alleles, const Allele& allele, VariantCaller::RefCallType refcall_type)
    {
        if (refcall_type == VariantCaller::RefCallType::Blocked && !alleles.empty() && are_adjacent(alleles.back(), allele)) {
            alleles.back() = Allele {
                get_encompassing(alleles.back(), allele),
                alleles.back().get_sequence() + allele.get_sequence()
            };
        } else {
            alleles.push_back(allele);
        }
    }
    
    // TODO: we should catch the case where an insertion has been called and push the refcall
    // block up a position, otherwise the returned reference allele (block) will never be called.
    std::vector<Allele>
    generate_candidate_reference_alleles(const std::vector<Allele>& callable_alleles,
                                         const std::vector<GenomicRegion>& called_regions,
                                         const std::vector<Variant>& candidates,
                                         VariantCaller::RefCallType refcall_type)
    {
        using std::cbegin; using std::cend; using std::advance;
        
        if (callable_alleles.empty() || refcall_type == VariantCaller::RefCallType::None) return {};
        
        if (candidates.empty()) return callable_alleles;
        
        auto allele_itr        = cbegin(callable_alleles);
        auto allele_end_itr    = cend(callable_alleles);
        auto called_itr        = cbegin(called_regions);
        auto called_end_itr    = cend(called_regions);
        auto candidate_itr     = cbegin(candidates);
        auto candidate_end_itr = cend(candidates);
        
        std::vector<Allele> result {};
        result.reserve(callable_alleles.size());
        
        while (allele_itr != allele_end_itr) {
            if (candidate_itr == candidate_end_itr) {
                append_allele(result, *allele_itr, refcall_type);
                std::copy(std::next(allele_itr), allele_end_itr, std::back_inserter(result));
                break;
            }
            
            if (called_itr == called_end_itr) {
                append_allele(result, *allele_itr, refcall_type);
                if (begins_before(*allele_itr, *candidate_itr)) {
                    advance(allele_itr, 1);
                } else {
                    advance(allele_itr, 2); // as next allele is alt from candidate
                    advance(candidate_itr, 1);
                }
            } else {
                if (is_same_region(*called_itr, *candidate_itr)) { // called candidate
                    while (is_before(*allele_itr, *called_itr)) {
                        append_allele(result, *allele_itr, refcall_type);
                        advance(allele_itr, 1);
                    }
                    advance(allele_itr, 2); // skip called variant
                    advance(candidate_itr, 1);
                    advance(called_itr, 1);
                } else if (begins_before(*called_itr, *candidate_itr)) { // parsimonised called candidate
                    if (!overlaps(*allele_itr, *called_itr)) {
                        append_allele(result, *allele_itr, refcall_type);
                        advance(allele_itr, 1);
                    } else {
                        if (begins_before(*allele_itr, *called_itr)) { // when variant has been left padded
                            append_allele(result, splice(*allele_itr, get_left_overhang(*allele_itr, *called_itr)), refcall_type);
                        }
                        
                        // skip contained alleles and candidates as they include called variants
                        allele_itr    = contained_range(allele_itr, allele_end_itr, *called_itr).end().base();
                        candidate_itr = contained_range(candidate_itr, candidate_end_itr, *called_itr).end().base();
                        
                        advance(called_itr, 1);
                    }
                } else {
                    append_allele(result, *allele_itr, refcall_type);
                    
                    if (begins_before(*allele_itr, *candidate_itr)) {
                        advance(allele_itr, 1);
                    } else {
                        advance(allele_itr, 2); // as next allele is alt from candidate
                        advance(candidate_itr, 1);
                    }
                }
            }
        }
        
        result.shrink_to_fit();
        
        return result;
    }
    
} // namespace Octopus
