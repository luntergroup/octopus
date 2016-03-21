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
#include <iostream>

#include "genomic_region.hpp"
#include "mappable.hpp"
#include "mappable_algorithms.hpp"
#include "mappable_flat_set.hpp"
#include "read_utils.hpp"
#include "haplotype.hpp"
#include "haplotype_generator.hpp"
#include "haplotype_liklihood_model.hpp"
#include "haplotype_filter.hpp"
#include "vcf_record.hpp"
#include "maths.hpp"
#include "logging.hpp"

#include "basic_haplotype_prior_model.hpp" // TODO: turn into factory

#include "timers.hpp"

namespace Octopus
{
// public methods
    
VariantCaller::CallerParameters::CallerParameters(unsigned max_haplotypes,
                                                  RefCallType refcall_type,
                                                  bool call_sites_only)
:
max_haplotypes {max_haplotypes},
refcall_type {refcall_type},
call_sites_only {call_sites_only}
{}

VariantCaller::VariantCaller(const ReferenceGenome& reference,
                             ReadPipe& read_pipe,
                             CandidateVariantGenerator&& candidate_generator,
                             std::unique_ptr<HaplotypePriorModel> haplotype_prior_model,
                             CallerParameters parameters)
:
reference_ {reference},
read_pipe_ {read_pipe},
refcall_type_ {parameters.refcall_type},
call_sites_only_ {parameters.call_sites_only},
max_haplotypes_ {parameters.max_haplotypes},
haplotype_prior_model_ {std::move(haplotype_prior_model)},
candidate_generator_ {std::move(candidate_generator)}
{}

auto generate_candidates(CandidateVariantGenerator& generator, const GenomicRegion& region,
                         const ReferenceGenome& reference)
{
    auto candidates = unique_left_align(generator.generate_candidates(region), reference);
    
    return MappableFlatSet<Variant> {
        std::make_move_iterator(std::begin(candidates)),
        std::make_move_iterator(std::end(candidates))
    };
}

VcfRecord annotate_record(VcfRecord::Builder& record, const ReadMap& reads)
{
    return record.build_once();
}

void append_annotated_calls(std::deque<VcfRecord>& curr_calls,
                            std::vector<VcfRecord::Builder>& new_calls,
                            const ReadMap& reads, const GenomicRegion& call_region)
{
    // Need to check calls are within call region as we may consider candidates
    // outside of this region when building haplotypes
    const auto it = std::find_if(std::begin(new_calls), std::end(new_calls),
                                 [&call_region] (const auto& call) {
                                     return region_begin(call_region) <= call.get_position() + 1;
                                 });
    
    const auto it2 = std::find_if(std::make_reverse_iterator(std::end(new_calls)),
                                  std::make_reverse_iterator(std::begin(new_calls)),
                                  [&call_region] (const auto& call) {
                                      return call.get_position() < region_end(call_region);
                                  }).base();
    
    //std::cout << "calling " << std::distance(it, it2) << " variants" << std::endl;
    
    std::transform(it, it2, std::back_inserter(curr_calls),
                   [&] (auto& record) { return annotate_record(record, reads); });
}

namespace debug
{
    template <typename S>
    void print_candidates(S&& stream, const MappableFlatSet<Variant>& candidates, bool number_only = false);
    void print_candidates(const MappableFlatSet<Variant>& candidates, bool number_only = false);
    
    template <typename S>
    void print_active_candidates(S&& stream, const MappableFlatSet<Variant>& candidates,
                                 const GenomicRegion& active_region, bool number_only = false);
    void print_active_candidates(const MappableFlatSet<Variant>& candidates,
                                 const GenomicRegion& active_region, bool number_only = false);
    
    template <typename S>
    void print_inactive_flanking_candidates(S&& stream, const MappableFlatSet<Variant>& candidates,
                                            const GenomicRegion& active_region,
                                            const GenomicRegion& haplotype_region,
                                            bool number_only = false);
    void print_inactive_flanking_candidates(const MappableFlatSet<Variant>& candidates,
                                            const GenomicRegion& active_region,
                                            const GenomicRegion& haplotype_region,
                                            bool number_only = false);
    
    enum class Resolution {Sequence, Alleles, VariantAlleles, SequenceAndAlleles, SequenceAndVariantAlleles};
    
    template <typename S>
    void print_haplotypes(S&& stream, const std::vector<Haplotype>& haplotypes,
                          Resolution resolution = Resolution::SequenceAndAlleles);
    void print_haplotypes(const std::vector<Haplotype>& haplotypes,
                          Resolution resolution = Resolution::SequenceAndAlleles);
    
    template <typename S, typename Map>
    void print_haplotype_posteriors(S&& stream, const Map& haplotype_posteriors, std::size_t n = 5);
    template <typename Map>
    void print_haplotype_posteriors(const Map& haplotype_posteriors, std::size_t n = 5);
}

auto copy_overlapped_to_vector(const MappableFlatSet<Variant>& candidates,
                               const GenomicRegion& region)
{
    const auto overlapped = overlap_range(candidates, region);
    return std::vector<Variant> {std::cbegin(overlapped), std::cend(overlapped)};
    //return copy_overlapped<std::vector<Variant>>(candidates, region);
}

const auto& haplotype_region(const std::vector<Haplotype>& haplotypes)
{
    return mapped_region(haplotypes.front());
}

void remove_passed_candidates(MappableFlatSet<Variant>& candidates,
                              const GenomicRegion& candidate_region,
                              const GenomicRegion& haplotype_region)
{
    if (begins_before(candidate_region, haplotype_region)) {
        const auto passed_region = left_overhang_region(candidate_region, haplotype_region);
        
        if (DEBUG_MODE) {
            Logging::DebugLogger log {};
            stream(log) << "Removing " << count_overlapped(candidates, passed_region)
                        << " passed candidates in region " << passed_region;
        }
        
        candidates.erase_overlapped(passed_region);
    }
}

//auto max_sequence_size(const std::vector<Haplotype>& haplotypes, const GenomicRegion& region)
//{
//    assert(!haplotypes.empty());
//    
//    std::vector<Haplotype::SizeType> region_sequence_sizes(haplotypes.size());
//    
//    std::transform(std::cbegin(haplotypes), std::cend(haplotypes),
//                   std::begin(region_sequence_sizes),
//                   [&region] (const auto& haplotype) {
//                       return haplotype.sequence_size(region);
//                   });
//    
//    return *std::max_element(std::cbegin(region_sequence_sizes), std::cend(region_sequence_sizes));
//}

auto calculate_flank_regions(const GenomicRegion& haplotype_region,
                             const GenomicRegion& active_region,
                             const MappableFlatSet<Variant>& candidates)
{
    auto lhs_flank = left_overhang_region(haplotype_region, active_region);
    auto rhs_flank = right_overhang_region(haplotype_region, active_region);
    
    const auto active_candidates = overlap_range(candidates, active_region);
    
    assert(!active_candidates.empty());
    
    if (is_empty_region(*leftmost_mappable(active_candidates)) && !is_empty_region(lhs_flank)) {
        lhs_flank = expand_rhs(lhs_flank, -1); // stops boundry insertions being considerd inactive
    }
    
    const auto lhs_inactive_candidates = overlap_range(candidates, lhs_flank);
    
    if (lhs_inactive_candidates.empty()) {
        lhs_flank = head_region(lhs_flank);
    } else {
        lhs_flank = closed_region(lhs_flank, rightmost_region(lhs_inactive_candidates));
    }
    
    if (is_empty_region(*rightmost_mappable(active_candidates)) && !is_empty_region(rhs_flank)) {
        rhs_flank = expand_lhs(rhs_flank, -1); // stops boundry insertions being considerd inactive
    }
    
    const auto rhs_inactive_candidates = overlap_range(candidates, rhs_flank);
    
    if (rhs_inactive_candidates.empty()) {
        rhs_flank = tail_region(rhs_flank);
    } else {
        rhs_flank = closed_region(leftmost_region(rhs_inactive_candidates), rhs_flank);
    }
    
    return std::make_pair(std::move(lhs_flank), std::move(rhs_flank));
}

auto calculate_flank_state(const std::vector<Haplotype>& haplotypes,
                           const GenomicRegion& active_region,
                           const MappableFlatSet<Variant>& candidates)
{
    const auto flanks = calculate_flank_regions(haplotype_region(haplotypes), active_region,
                                                candidates);
    
    return HaplotypeLikelihoodModel::FlankState {
        active_region.get_contig_region(),
        flanks.first.get_contig_region(),
        flanks.second.get_contig_region()
    };
}

bool all_empty(const ReadMap& reads)
{
    return std::all_of(std::cbegin(reads), std::cend(reads),
                       [] (const auto& p) { return p.second.empty(); });
}

auto calculate_candidate_region(const GenomicRegion& call_region, const ReadMap& reads,
                                const CandidateVariantGenerator& candidate_generator)
{
    if (!candidate_generator.requires_reads()) return call_region;
    
    if (all_empty(reads)) {
        return call_region;
    }
    
    return encompassing_region(reads);
}

bool have_passed(const GenomicRegion& next_active_region, const GenomicRegion& active_region)
{
    return is_after(next_active_region, active_region) && active_region != next_active_region;
}

bool has_moved_forward(const GenomicRegion& next_active_region, const GenomicRegion& active_region)
{
    return begins_before(active_region, next_active_region) || active_region == next_active_region;
}

std::deque<VcfRecord> VariantCaller::call_variants(const GenomicRegion& call_region,
                                                   ProgressMeter& progress_meter) const
{
    Logging::DebugLogger debug_log {};
    
    resume_timer(init_timer);
    
    ReadMap reads;
    
    std::deque<VcfRecord> result {};
    
    if (candidate_generator_.requires_reads()) {
        reads = read_pipe_.get().fetch_reads(call_region);
        
        add_reads(reads, candidate_generator_);
        
        if (!refcalls_requested() && all_empty(reads)) {
            return result;
        }
        
        if (DEBUG_MODE) {
            stream(debug_log) << "There are " << count_reads(reads) << " reads";
        }
    }
    
    const auto candidate_region = calculate_candidate_region(call_region, reads, candidate_generator_);
    
    if (DEBUG_MODE) {
        stream(debug_log) << "Generating candidates in region " << candidate_region;
    }
    
    auto candidates = generate_candidates(candidate_generator_, candidate_region, reference_);
    
    candidate_generator_.clear();
    
    if (DEBUG_MODE) {
        debug::print_candidates(stream(debug_log), candidates);
    }
    
    if (!refcalls_requested() && candidates.empty()) {
        return result;
    }
    
    if (!candidate_generator_.requires_reads()) {
        // as we didn't fetch them earlier
        reads = read_pipe_.get().fetch_reads(extract_regions(candidates));
    }
    
    constexpr unsigned max_indicators {0};
    constexpr double min_phase_score {0.99};
    
    HaplotypeGenerator generator {candidate_region, reference_, candidates, reads, max_haplotypes_, max_indicators};
    Phaser phaser {min_phase_score};
    
    std::vector<Haplotype> haplotypes;
    GenomicRegion active_region;
    
    auto completed_region = head_region(call_region);
    
    pause_timer(init_timer);
    
    const auto& samples = read_pipe_.get().get_samples();
    
    HaplotypeLikelihoodCache haplotype_likelihoods {max_haplotypes_, samples};
    
    while (true) {
        resume_timer(haplotype_generation_timer);
        std::tie(haplotypes, active_region) = generator.progress();
        pause_timer(haplotype_generation_timer);
        
        if (is_after(active_region, call_region) || haplotypes.empty()) {
            progress_meter.log_completed(active_region);
            break;
        }
        
        remove_passed_candidates(candidates, candidate_region, haplotype_region(haplotypes));
        
        if (DEBUG_MODE) {
            stream(debug_log) << "Active region is " << active_region;
            debug::print_active_candidates(stream(debug_log), candidates, active_region);
            stream(debug_log) << "Haplotype region is " << haplotype_region(haplotypes);
            debug::print_inactive_flanking_candidates(stream(debug_log), candidates, active_region,
                                                      haplotype_region(haplotypes));
        }
        
        const auto active_reads = copy_overlapped(reads, active_region);
        
        if (DEBUG_MODE) {
            stream(debug_log) << "There are " << haplotypes.size() << " initial haplotypes";
            stream(debug_log) << "There are " << count_reads(active_reads) << " active reads";
        }
        
        resume_timer(likelihood_timer);
        haplotype_likelihoods.populate(active_reads, haplotypes,
                                        calculate_flank_state(haplotypes, active_region, candidates));
        pause_timer(likelihood_timer);
        
        resume_timer(haplotype_fitler_timer);
        auto removed_haplotypes = filter_n_haplotypes(haplotypes, samples,
                                                      haplotype_likelihoods, max_haplotypes_);
        pause_timer(haplotype_fitler_timer);
        
        if (haplotypes.empty()) {
            if (DEBUG_MODE) debug_log << "Filtered all haplotypes";
            // This can only happen if all haplotypes have equal likelihood.
            // TODO: is there anything else we can do?
            generator.remove_haplotypes(removed_haplotypes);
            continue;
        }
        
        if (TRACE_MODE) {
            Logging::TraceLogger trace_log {};
            stream(trace_log) << "Filtered " << removed_haplotypes.size() << " haplotypes:";
            debug::print_haplotypes(stream(trace_log), removed_haplotypes);
        } else if (DEBUG_MODE) {
            stream(debug_log) << "Filtered " << removed_haplotypes.size() << " haplotypes";
        }
        
        resume_timer(likelihood_timer);
        haplotype_likelihoods.erase(removed_haplotypes);
        pause_timer(likelihood_timer);
        
        // Compute haplotype priors after likelihood filtering as prior model may use
        // interdependencies between haplotypes.
        resume_timer(prior_model_timer);
        auto haplotype_priors = haplotype_prior_model_->compute_maximum_entropy_haplotype_set(haplotypes);
        pause_timer(prior_model_timer);
        
        if (TRACE_MODE) {
            Logging::TraceLogger trace_log {};
            debug::print_haplotype_priors(stream(trace_log), haplotype_priors, -1);
        } else if (DEBUG_MODE) {
            debug::print_haplotype_priors(stream(debug_log), haplotype_priors);
        }
        
        resume_timer(haplotype_generation_timer);
        generator.keep_haplotypes(haplotypes);
        generator.remove_haplotypes(removed_haplotypes);
        removed_haplotypes.clear();
        pause_timer(haplotype_generation_timer);
        
        if (TRACE_MODE) {
            Logging::TraceLogger trace_log {};
            stream(trace_log) << "There are " << haplotypes.size() << " final haplotypes";
            debug::print_read_haplotype_liklihoods(stream(trace_log), haplotypes, active_reads,
                                                   haplotype_likelihoods, -1);
        }
        
        resume_timer(latent_timer);
        const auto caller_latents = infer_latents(samples, haplotypes, haplotype_priors,
                                                  haplotype_likelihoods);
        pause_timer(latent_timer);
        
        haplotype_priors.clear();
        haplotype_likelihoods.clear();
        
        if (TRACE_MODE) {
            Logging::TraceLogger trace_log {};
            debug::print_haplotype_posteriors(stream(trace_log),
                                              *caller_latents->get_haplotype_posteriors(), -1);
        } else if (DEBUG_MODE) {
            debug::print_haplotype_posteriors(stream(debug_log),
                                              *caller_latents->get_haplotype_posteriors());
        }
        
        resume_timer(phasing_timer);
        const auto phase_set = phaser.try_phase(haplotypes, *caller_latents->get_genotype_posteriors(),
                                                copy_overlapped_to_vector(candidates, haplotype_region(haplotypes)));
        pause_timer(phasing_timer);
        
        if (DEBUG_MODE) {
            if (phase_set) {
                debug::print_phase_sets(stream(debug_log), *phase_set);
            } else {
                debug_log << "No partial phasings found";
            }
        }
        
        auto unphased_active_region = active_region;
        
        if (overlaps(active_region, call_region) && phase_set) {
            assert(!is_empty_region(phase_set->region));
            
            if (DEBUG_MODE) stream(debug_log) << "Phased region is " << phase_set->region;
            
            auto active_candidates = copy_overlapped_to_vector(candidates, phase_set->region);
            
            resume_timer(allele_generator_timer);
            auto alleles = generate_callable_alleles(phase_set->region, active_candidates,
                                                     refcall_type_, reference_);
            pause_timer(allele_generator_timer);
            
            resume_timer(calling_timer);
            auto curr_results = call_variants(active_candidates, alleles, caller_latents.get(),
                                              *phase_set, active_reads);
            
            append_annotated_calls(result, curr_results, active_reads, call_region);
            pause_timer(calling_timer);
            
            auto remaining_active_region = right_overhang_region(active_region, phase_set->region);
            
            resume_timer(haplotype_generation_timer);
            generator.force_forward(remaining_active_region);
            pause_timer(haplotype_generation_timer);
            
            unphased_active_region = right_overhang_region(active_region, phase_set->region);
        }
        
        resume_timer(haplotype_generation_timer);
        auto next_active_region = generator.tell_next_active_region();
        pause_timer(haplotype_generation_timer);
        
        if (!have_passed(next_active_region, active_region)) {
            auto removable_haplotypes = get_removable_haplotypes(haplotypes,
                                                                 *caller_latents->get_haplotype_posteriors(),
                                                                 unphased_active_region);
            
            if (TRACE_MODE) {
                Logging::TraceLogger trace_log {};
                stream(trace_log) << "Removing " << removable_haplotypes.size() << " haplotypes:";
                debug::print_haplotypes(stream(trace_log), removable_haplotypes);
            } else if (DEBUG_MODE) {
                stream(debug_log) << "Removing " << removable_haplotypes.size() << " haplotypes";
            }
            
            resume_timer(haplotype_generation_timer);
            generator.remove_haplotypes(removable_haplotypes);
            next_active_region = generator.tell_next_active_region();
            pause_timer(haplotype_generation_timer);
        }
        
        if (overlaps(active_region, call_region) && has_moved_forward(next_active_region, active_region)) {
            auto passed_region   = left_overhang_region(active_region, next_active_region);
            auto uncalled_region = overlapped_region(active_region, passed_region);
            
            if (phase_set && ends_before(phase_set->region, passed_region)) {
                uncalled_region = right_overhang_region(passed_region, phase_set->region);
            }
            
            const auto active_candidates = copy_overlapped_to_vector(candidates, uncalled_region);
            
            resume_timer(phasing_timer);
            
            if (DEBUG_MODE) {
                stream(debug_log) << "Force phasing " << uncalled_region;
            }
            
            const auto phasings = phaser.force_phase(haplotypes,
                                                     *caller_latents->get_genotype_posteriors(),
                                                     active_candidates);
            pause_timer(phasing_timer);
            
            if (DEBUG_MODE) {
                debug::print_phase_sets(stream(debug_log), phasings);
            }
            
            haplotypes.clear();
            haplotypes.shrink_to_fit();
            
            resume_timer(allele_generator_timer);
            //uncalled_region = call_region; // DEBUG
            auto alleles = generate_callable_alleles(uncalled_region, active_candidates,
                                                     refcall_type_, reference_);
            pause_timer(allele_generator_timer);
            
            resume_timer(calling_timer);
            auto curr_results = call_variants(active_candidates, alleles, caller_latents.get(),
                                              phasings, active_reads);
            
            append_annotated_calls(result, curr_results, active_reads, call_region);
            pause_timer(calling_timer);
            
            completed_region = encompassing_region(completed_region, passed_region);
        }
        
        progress_meter.log_completed(active_region);
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
        if (p.second < 1e-12) {
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
        std::merge(make_move_iterator(begin(reference_alleles)),
                   make_move_iterator(end(reference_alleles)),
                   make_move_iterator(begin(variant_alleles)),
                   make_move_iterator(end(variant_alleles)),
                   back_inserter(result));
    } else {
        result.reserve(variant_alleles.size() + sum_region_sizes(uncovered_regions));
        
        auto uncovered_itr = begin(uncovered_regions);
        auto uncovered_end = end(uncovered_regions);
        
        for (auto&& variant_allele : variant_alleles) {
            if (uncovered_itr != uncovered_end && begins_before(*uncovered_itr, variant_allele)) {
                auto alleles = make_positional_reference_alleles(*uncovered_itr, reference);
                result.insert(end(result),
                              make_move_iterator(begin(alleles)),
                              make_move_iterator(end(alleles)));
                std::advance(uncovered_itr, 1);
            }
            result.push_back(std::move(variant_allele));
        }
        
        if (uncovered_itr != uncovered_end) {
            auto alleles = make_positional_reference_alleles(*uncovered_itr, reference);
            result.insert(end(result),
                          make_move_iterator(begin(alleles)),
                          make_move_iterator(end(alleles)));
        }
    }
    
    return result;
}

template <typename ForwardIt>
ForwardIt find_next(ForwardIt first, ForwardIt last, const Variant& candidate)
{
    return std::find_if_not(first, last,
                            [&] (const Allele& allele) {
                                return is_same_region(allele, candidate);
                            });
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
    using std::cbegin; using std::cend;
    
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
                ++allele_itr;
            } else {
                allele_itr = find_next(allele_itr, allele_end_itr, *candidate_itr);
                ++candidate_itr;
            }
        } else {
            if (is_same_region(*called_itr, *candidate_itr)) { // called candidate
                while (is_before(*allele_itr, *called_itr)) {
                    append_allele(result, *allele_itr, refcall_type);
                    ++allele_itr;
                }
                allele_itr = find_next(allele_itr, allele_end_itr, *candidate_itr);
                ++candidate_itr;
                ++called_itr;
            } else if (begins_before(*called_itr, *candidate_itr)) { // parsimonised called candidate
                if (!overlaps(*allele_itr, *called_itr)) {
                    append_allele(result, *allele_itr, refcall_type);
                    ++allele_itr;
                } else {
                    if (begins_before(*allele_itr, *called_itr)) { // when variant has been left padded
                        append_allele(result, splice(*allele_itr, left_overhang_region(*allele_itr, *called_itr)),
                                      refcall_type);
                    }
                    
                    // skip contained alleles and candidates as they include called variants
                    allele_itr    = cend(contained_range(allele_itr, allele_end_itr, *called_itr)).base();
                    candidate_itr = cend(contained_range(candidate_itr, candidate_end_itr, *called_itr)).base();
                    
                    ++called_itr;
                }
            } else {
                append_allele(result, *allele_itr, refcall_type);
                
                if (begins_before(*allele_itr, *candidate_itr)) {
                    ++allele_itr;
                } else {
                    allele_itr = find_next(allele_itr, allele_end_itr, *candidate_itr);
                    ++candidate_itr;
                }
            }
        }
    }
    
    result.shrink_to_fit();
    
    return result;
}

namespace debug
{
    template <typename S>
    void print_candidates(S&& stream, const MappableFlatSet<Variant>& candidates, bool number_only)
    {
        if (candidates.empty()) {
            stream << "There are no candidates" << '\n';
        } else {
            if (candidates.size() == 1) {
                stream << "There is 1 candidate" << '\n';
            } else {
                stream << "There are " << candidates.size() << " candidates: " << '\n';
            }
            
            if (!number_only) {
                for (const auto& c : candidates) stream << c << '\n';
            }
        }
    }
    
    void print_candidates(const MappableFlatSet<Variant>& candidates, bool number_only)
    {
        print_candidates(std::cout, candidates, number_only);
    }
    
    template <typename S>
    void print_active_candidates(S&& stream, const MappableFlatSet<Variant>& candidates,
                                 const GenomicRegion& active_region, bool number_only)
    {
        const auto active_candidates = overlap_range(candidates, active_region);
        
        if (active_candidates.empty()) {
            stream << "There are no active candidates" << '\n';
        } else {
            const auto num_active = size(active_candidates);
            
            if (num_active == 1) {
                stream << "There is 1 active candidate" << '\n';
            } else {
                stream << "There are " << num_active << " active candidates: " << '\n';
            }
            
            if (!number_only) {
                for (const auto& c : active_candidates) stream << c << '\n';
            }
        }
    }
    
    void print_active_candidates(const MappableFlatSet<Variant>& candidates,
                                 const GenomicRegion& active_region, bool number_only)
    {
        print_active_candidates(std::cout, candidates, active_region, number_only);
    }
    
    template <typename S>
    void print_inactive_flanking_candidates(S&& stream, const MappableFlatSet<Variant>& candidates,
                                            const GenomicRegion& active_region,
                                            const GenomicRegion& haplotype_region,
                                            bool number_only)
    {
        const auto flanks = calculate_flank_regions(haplotype_region, active_region, candidates);
        
        stream << "Haplotype flank regions are " << flanks.first << " " << flanks.second << '\n';
        
        const auto lhs_inactive_candidates = overlap_range(candidates, flanks.first);
        const auto rhs_inactive_candidates = overlap_range(candidates, flanks.second);
        
        const auto num_lhs_inactives = size(lhs_inactive_candidates);
        const auto num_rhs_inactives = size(rhs_inactive_candidates);
        
        if (lhs_inactive_candidates.empty()) {
            if (rhs_inactive_candidates.empty()) {
                stream << "There are no inactive flanking candidates" << '\n';
                return;
            }
            
            stream << "There are no lhs inactive flanking candidates" << '\n';
            
            if (num_rhs_inactives == 1) {
                stream << "There is 1 rhs inactive flanking candidates: " << '\n';
            } else {
                stream << "There are " << num_rhs_inactives << " rhs inactive flanking candidates: " << '\n';
            }
            
            if (!number_only) {
                for (const auto& c : rhs_inactive_candidates) stream << c << '\n';
            }
            
            return;
        }
        
        if (num_lhs_inactives == 1) {
            stream << "There is 1 lhs inactive flanking candidates: " << '\n';
        } else {
            stream << "There are " << num_lhs_inactives << " lhs inactive flanking candidates: " << '\n';
        }
        
        if (!number_only) {
            for (const auto& c : lhs_inactive_candidates) stream << c << '\n';
        }
        
        if (rhs_inactive_candidates.empty()) {
            stream << "There are no rhs inactive flanking candidates" << '\n';
        } else {
            if (num_rhs_inactives == 1) {
                stream << "There is 1 rhs inactive flanking candidates: " << '\n';
            } else {
                stream << "There are " << num_rhs_inactives << " rhs inactive flanking candidates: " << '\n';
            }
            
            if (!number_only) {
                for (const auto& c : rhs_inactive_candidates) stream << c << '\n';
            }
        }
    }
    
    void print_inactive_flanking_candidates(const MappableFlatSet<Variant>& candidates,
                                            const GenomicRegion& active_region,
                                            const GenomicRegion& haplotype_region,
                                            bool number_only)
    {
        print_inactive_flanking_candidates(std::cout, candidates, active_region, haplotype_region, number_only);
    }
    
    template <typename S>
    void print_haplotypes(S&& stream, const std::vector<Haplotype>& haplotypes,
                          Resolution resolution)
    {
        stream << "Printing " << haplotypes.size() << " haplotypes" << '\n';
        
        for (const auto& haplotype : haplotypes) {
            if (resolution == Resolution::Sequence || resolution == Resolution::SequenceAndAlleles
                || resolution == Resolution::SequenceAndVariantAlleles) {
                stream << haplotype << '\n';
            }
            if (resolution == Resolution::Alleles || resolution == Resolution::SequenceAndAlleles) {
                ::debug::print_alleles(haplotype); stream << '\n';
            } else if (resolution != Resolution::Sequence) {
                ::debug::print_variant_alleles(haplotype); stream << '\n';
            }
        }
    }
    
    void print_haplotypes(const std::vector<Haplotype>& haplotypes, const Resolution resolution)
    {
        print_haplotypes(std::cout, haplotypes, resolution);
    }
    
    template <typename S, typename Map>
    void print_haplotype_posteriors(S&& stream, const Map& haplotype_posteriors, std::size_t n)
    {
        const auto m = std::min(haplotype_posteriors.size(), n);
        
        if (m == haplotype_posteriors.size()) {
            stream << "Printing top " << m << " haplotype posteriors" << '\n';
        } else {
            stream << "Printing all haplotype posteriors" << '\n';
        }
        
        std::vector<std::pair<std::reference_wrapper<const Haplotype>, double>> v {};
        v.reserve(haplotype_posteriors.size());
        
        std::copy(std::cbegin(haplotype_posteriors), std::cend(haplotype_posteriors),
                  std::back_inserter(v));
        
        const auto mth = std::next(std::begin(v), m);
        
        std::partial_sort(std::begin(v), mth, std::end(v),
                          [] (const auto& lhs, const auto& rhs) {
                              return lhs.second > rhs.second;
                          });
        
        std::for_each(std::begin(v), mth,
                      [&] (const auto& p) {
                          ::debug::print_variant_alleles(stream, p.first);
                          stream << " " << p.second << '\n';
                      });
    }
    
    template <typename Map>
    void print_haplotype_posteriors(const Map& haplotype_posteriors, std::size_t n)
    {
        print_haplotype_posteriors(haplotype_posteriors, n);
    }
}
} // namespace Octopus
