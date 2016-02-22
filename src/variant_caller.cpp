//
//  variant_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 15/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant_caller.hpp"

#include <algorithm>
#include <utility>
#include <tuple>
#include <iterator>
#include <cassert>

#include "genomic_region.hpp"
#include "mappable.hpp"
#include "mappable_algorithms.hpp"
#include "read_utils.hpp"
#include "mappable_algorithms.hpp"
#include "haplotype.hpp"
#include "vcf_record.hpp"

#include "basic_haplotype_prior_model.hpp" // TODO: turn into factory

#include "haplotype_generator.hpp"

#include <iostream> // TEST

namespace Octopus
{
// public methods

VariantCaller::VariantCaller(const ReferenceGenome& reference,
                             ReadPipe& read_pipe,
                             CandidateVariantGenerator&& candidate_generator,
                             RefCallType refcall_type)
:
reference_ {reference},
read_pipe_ {read_pipe},
refcall_type_ {refcall_type},
haplotype_prior_model_ {std::make_unique<BasicHaplotypePriorModel>(reference_)},
candidate_generator_ {std::move(candidate_generator)}
{}

VariantCaller::VariantCaller(const ReferenceGenome& reference,
                             ReadPipe& read_pipe,
                             CandidateVariantGenerator&& candidate_generator,
                             std::unique_ptr<HaplotypePriorModel> haplotype_prior_model,
                             RefCallType refcall_type)
:
reference_ {reference},
read_pipe_ {read_pipe},
refcall_type_ {refcall_type},
haplotype_prior_model_ {std::move(haplotype_prior_model)},
candidate_generator_ {std::move(candidate_generator)}
{}

size_t VariantCaller::num_buffered_reads() const noexcept
{
    return 0;
}

VcfRecord annotate_record(VcfRecord::Builder& record, const ReadMap& reads)
{
    return record.build_once();
}

void append_annotated_calls(std::vector<VcfRecord>& curr_calls,
                            std::vector<VcfRecord::Builder>& new_calls,
                            const ReadMap& reads)
{
    curr_calls.reserve(curr_calls.size() + new_calls.size());
    std::transform(std::begin(new_calls), std::end(new_calls), std::back_inserter(curr_calls),
                   [&] (auto& record) { return annotate_record(record, reads); });
}

namespace debug
{
    void print_candidates(const std::vector<Variant>& candidates);
    enum class Resolution {Sequence, Alleles, VariantAlleles, SequenceAndAlleles, SequenceAndVariantAlleles};
    void print_haplotypes(const std::vector<Haplotype>& haplotypes,
                          Resolution resolution = Resolution::SequenceAndAlleles);
}

std::vector<Haplotype>
filter_haplotypes(std::vector<Haplotype>& haplotypes, const ReadMap& reads,
                  const HaplotypeLikelihoodCache& haplotype_likelihoods, const size_t n)
{
    std::vector<Haplotype> result {};
    
    if (haplotypes.size() <= n) return result;
    
    std::unordered_map<Haplotype, double> max_liklihoods {haplotypes.size()};
    
    for (const auto& haplotype : haplotypes) {
        double max_read_liklihood {0};
        
        for (const auto& sample_reads : reads) {
            for (const auto& read : sample_reads.second) {
                const auto cur_read_liklihood = haplotype_likelihoods.log_probability(read, haplotype);
                if (cur_read_liklihood > max_read_liklihood) max_read_liklihood = cur_read_liklihood;
            }
        }
        
        max_liklihoods.emplace(haplotype, max_read_liklihood);
    }
    
    const auto nth = std::next(std::begin(haplotypes), n);
    
    std::nth_element(std::begin(haplotypes), nth, std::end(haplotypes),
                     [&] (const auto& lhs, const auto& rhs) {
                         return max_liklihoods.at(lhs) > max_liklihoods.at(rhs);
                     });
    
    result.assign(nth, std::end(haplotypes));
    
    haplotypes.erase(nth, std::end(haplotypes));
    
    return result;
}

std::vector<VcfRecord> VariantCaller::call_variants(const GenomicRegion& call_region) const
{
    assert(!is_empty_region(call_region));
    
    ReadMap reads;
    
    if (candidate_generator_.requires_reads()) {
        reads = read_pipe_.get().fetch_reads(call_region);
        add_reads(reads, candidate_generator_);
    }
    
    const auto candidates = unique_left_align(candidate_generator_.get_candidates(call_region), reference_);
    
    candidate_generator_.clear();
    
    debug::print_candidates(candidates);
    
    std::vector<VcfRecord> result {};
    
    if (candidates.empty() && !refcalls_requested()) {
        return result;
    }
    
    if (!candidate_generator_.requires_reads()) {
        // TODO: we could be more selective and only fetch reads overlapping candidates
        reads = read_pipe_.get().fetch_reads(call_region);
    }
    
    constexpr unsigned max_haplotypes {64};
    constexpr unsigned max_indicators {0};
    constexpr double min_phase_score {0.99};
    
    HaplotypeGenerator generator {call_region, reference_, candidates, reads, max_haplotypes, max_indicators};
    Phaser phaser {min_phase_score};
    
    std::vector<Haplotype> haplotypes;
    GenomicRegion active_region;
    
    while (true) {
        std::tie(haplotypes, active_region) = generator.progress();
        
        if (haplotypes.empty()) break;
        
        std::cout << "active region is " << active_region << '\n';
        
        const auto active_region_reads = copy_overlapped(reads, active_region);
        
        std::cout << "there are " << count_reads(active_region_reads) << " reads in active region" << '\n';
        std::cout << "active_region_reads region is " << encompassing_region(active_region_reads) << '\n';
        
        std::cout << "there are " << haplotypes.size() << " initial haplotypes" << '\n';
        
        HaplotypeLikelihoodCache haplotype_likelihoods {active_region_reads, haplotypes};
        
        std::cout << "computed haplotype likelihoods" << '\n';
        
        auto removed_haplotypes = filter_haplotypes(haplotypes, active_region_reads, haplotype_likelihoods,
                                                    max_haplotypes);
        
        // Compute haplotype priors after likelihood filtering as prior model may use
        // interdependencies between haplotypes.
        auto haplotype_priors = haplotype_prior_model_->compute_maximum_entropy_haplotype_set(haplotypes);
        
        generator.keep_haplotypes(haplotypes);
        generator.remove_haplotypes(removed_haplotypes);
        
        std::cout << "there are " << haplotypes.size() << " final haplotypes" << '\n';
        
        const auto caller_latents = infer_latents(haplotypes, haplotype_priors,
                                                  haplotype_likelihoods, active_region_reads);
        
        haplotype_priors.clear();
        haplotype_likelihoods.clear();
        
        const auto phase_set = phaser.try_phase(haplotypes, *caller_latents->get_genotype_posteriors());
        
        auto unphased_active_region = active_region;
        
        if (phase_set) {
            assert(!is_empty_region(phase_set->region));
            
            std::cout << "phased region is " << phase_set->region << '\n';
            
            auto active_candidates = copy_overlapped(candidates, phase_set->region);
            
            auto alleles = generate_callable_alleles(phase_set->region, active_candidates,
                                                     refcall_type_, reference_);
            
            auto curr_results = call_variants(active_candidates, alleles, caller_latents.get(),
                                              *phase_set, active_region_reads);
            
            append_annotated_calls(result, curr_results, active_region_reads);
            
            auto remaining_active_region = right_overhang_region(active_region, phase_set->region);
            
            generator.force_forward(remaining_active_region);
            
            unphased_active_region = right_overhang_region(active_region, phase_set->region);
        }
        
        // voluntarily remove some haplotypes in unphased_active_region
        auto removable_haplotypes = get_removable_haplotypes(haplotypes, *caller_latents->get_haplotype_posteriors(),
                                                             unphased_active_region);
        generator.remove_haplotypes(removable_haplotypes);
        
        auto next_active_region = generator.tell_next_active_region();
        
        if (begins_before(active_region, next_active_region)) {
            auto passed_region = left_overhang_region(active_region, next_active_region);
            
            auto uncalled_region = passed_region;
            
            if (phase_set && ends_before(phase_set->region, passed_region)) {
                uncalled_region = right_overhang_region(passed_region, phase_set->region);
            }
            
            if (!is_empty_region(uncalled_region)) {
                const auto forced_phasing = phaser.force_phase(haplotypes, *caller_latents->get_genotype_posteriors());
                
                auto active_candidates = copy_overlapped(candidates, uncalled_region);
                
                auto alleles = generate_callable_alleles(uncalled_region, active_candidates,
                                                         refcall_type_, reference_);
                
                auto curr_results = call_variants(active_candidates, alleles, caller_latents.get(),
                                                  forced_phasing, active_region_reads);
                
                append_annotated_calls(result, curr_results, active_region_reads);
            }
        }
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
    return is_empty_region(region);
}

std::vector<Haplotype>
VariantCaller::get_removable_haplotypes(const std::vector<Haplotype>& haplotypes,
                         const CallerLatents::HaplotypePosteiorMap& haplotype_posteriors,
                         const GenomicRegion& region) const
{
    assert(!haplotypes.empty() && contains(haplotypes.front().get_region(), region));
    
    std::vector<Haplotype> result {};
    
    // TODO
    for (const auto& p : haplotype_posteriors) {
        if (p.second < 1e-5) {
            result.push_back(p.first);
        }
    }
    
    return result;
}

// non-member methods

std::vector<Allele>
generate_callable_alleles(const GenomicRegion& region,
                          const std::vector<Variant>& variants,
                          const VariantCaller::RefCallType refcall_type,
                          const ReferenceGenome& reference)
{
    using std::begin; using std::end; using std::make_move_iterator; using std::back_inserter;
    
    auto overlapped_variants = copy_overlapped(variants, region);
    
    if (is_empty_region(region) && overlapped_variants.empty()) return std::vector<Allele> {};
    
    if (overlapped_variants.empty()) {
        switch (refcall_type) {
            case VariantCaller::RefCallType::Positional:
                return make_positional_reference_alleles(region, reference);
            case VariantCaller::RefCallType::Blocked:
                return std::vector<Allele> {make_reference_allele(region, reference)};
            default:
                return {};
        }
    }
    
    auto variant_alleles = decompose(overlapped_variants);
    
    if (refcall_type == VariantCaller::RefCallType::None) return variant_alleles;
    
    auto covered_regions   = extract_covered_regions(overlapped_variants);
    auto uncovered_regions = extract_intervening_regions(covered_regions, region);
    
    std::vector<Allele> result {};
    
    if (refcall_type == VariantCaller::RefCallType::Blocked) {
        auto reference_alleles = make_reference_alleles(uncovered_regions, reference);
        result.reserve(reference_alleles.size() + variant_alleles.size());
        std::merge(make_move_iterator(begin(reference_alleles)),make_move_iterator(end(reference_alleles)),
                   make_move_iterator(begin(variant_alleles)), make_move_iterator(end(variant_alleles)),
                   back_inserter(result));
    } else {
        result.reserve(variant_alleles.size() + sum_region_sizes(uncovered_regions));
        
        auto uncovered_itr = begin(uncovered_regions);
        auto uncovered_end = end(uncovered_regions);
        
        for (auto&& variant_allele : variant_alleles) {
            if (uncovered_itr != uncovered_end && begins_before(*uncovered_itr, variant_allele)) {
                auto alleles = make_positional_reference_alleles(*uncovered_itr, reference);
                result.insert(end(result), make_move_iterator(begin(alleles)),
                              make_move_iterator(end(alleles)));
                std::advance(uncovered_itr, 1);
            }
            result.push_back(std::move(variant_allele));
        }
        
        if (uncovered_itr != uncovered_end) {
            auto alleles = make_positional_reference_alleles(*uncovered_itr, reference);
            result.insert(end(result), make_move_iterator(begin(alleles)),
                          make_move_iterator(end(alleles)));
        }
    }
    
    return result;
}

void append_allele(std::vector<Allele>& alleles, const Allele& allele,
                   const VariantCaller::RefCallType refcall_type)
{
    if (refcall_type == VariantCaller::RefCallType::Blocked && !alleles.empty()
        && are_adjacent(alleles.back(), allele)) {
        alleles.back() = Allele {encompassing_region(alleles.back(), allele),
                                    alleles.back().get_sequence() + allele.get_sequence()};
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
                                     const VariantCaller::RefCallType refcall_type)
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
                        append_allele(result, splice(*allele_itr, left_overhang_region(*allele_itr, *called_itr)), refcall_type);
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
    
namespace debug
{
    void print_candidates(const std::vector<Variant>& candidates)
    {
        if (candidates.empty()) {
            std::cout << "there are no candidates" << '\n';
        } else {
            std::cout << "found " << candidates.size() << " candidates: " << '\n';
            for (const auto& c : candidates) std::cout << c << '\n';
        }
    }
    
    void print_haplotypes(const std::vector<Haplotype>& haplotypes, const Resolution resolution)
    {
        std::cout << "printing " << haplotypes.size() << " haplotypes" << '\n';
        for (const auto& haplotype : haplotypes) {
            if (resolution == Resolution::Sequence || resolution == Resolution::SequenceAndAlleles
                || resolution == Resolution::SequenceAndVariantAlleles) {
                std::cout << haplotype << '\n';
            }
            if (resolution == Resolution::Alleles || resolution == Resolution::SequenceAndAlleles) {
                print_alleles(haplotype); std::cout << '\n';
            } else if (resolution != Resolution::Sequence) {
                print_variant_alleles(haplotype); std::cout << '\n';
            }
        }
    }
}
} // namespace Octopus
