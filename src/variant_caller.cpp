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
#include "haplotype_generator.hpp"
#include "haplotype_liklihood_model.hpp"
#include "haplotype_filter.hpp"
#include "vcf_record.hpp"
#include "maths.hpp"

#include "basic_haplotype_prior_model.hpp" // TODO: turn into factory

#include "logging.hpp"

#include <iostream> // DEBUG
#include "timers.hpp"

namespace Octopus
{
// public methods

VariantCaller::VariantCaller(const ReferenceGenome& reference,
                             ReadPipe& read_pipe,
                             CandidateVariantGenerator&& candidate_generator,
                             const unsigned max_haplotypes,
                             std::unique_ptr<HaplotypePriorModel> haplotype_prior_model,
                             RefCallType refcall_type)
:
reference_ {reference},
read_pipe_ {read_pipe},
refcall_type_ {refcall_type},
max_haplotypes_ {max_haplotypes},
haplotype_prior_model_ {std::move(haplotype_prior_model)},
candidate_generator_ {std::move(candidate_generator)}
{}

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
    void print_candidates(const std::vector<Variant>& candidates, bool number_only = false);
    enum class Resolution {Sequence, Alleles, VariantAlleles, SequenceAndAlleles, SequenceAndVariantAlleles};
    void print_haplotypes(const std::vector<Haplotype>& haplotypes,
                          Resolution resolution = Resolution::SequenceAndAlleles);
    void print_progress(const GenomicRegion& last_active_region,
                        const GenomicRegion& call_region,
                        GenomicRegion::SizeType step_size);
    
    template <typename Map>
    void print_haplotype_posteriors(const Map& haplotype_posteriors, std::size_t n = 20);
}

auto flank_state(const std::vector<Haplotype>& haplotypes,
                 const GenomicRegion& active_region,
                 const std::vector<Variant>& candidates)
{
    const auto& haplotype_region = haplotypes.front().get_region();
    
    return HaplotypeLikelihoodModel::FlankState {
        active_region.get_contig_region(),
        has_overlapped(candidates, left_overhang_region(haplotype_region, active_region)),
        has_overlapped(candidates, right_overhang_region(haplotype_region, active_region))
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

std::deque<VcfRecord> VariantCaller::call_variants(const GenomicRegion& call_region) const
{
    resume_timer(init_timer);
    
    ReadMap reads;
    
    std::deque<VcfRecord> result {};
    
    if (candidate_generator_.requires_reads()) {
        reads = read_pipe_.get().fetch_reads(call_region);
        
        add_reads(reads, candidate_generator_);
        
        if (!refcalls_requested() && all_empty(reads)) {
            return result;
        }
        
        //std::cout << "there are " << count_reads(reads) << " reads" << '\n';
    }
    
    const auto candidate_region = calculate_candidate_region(call_region, reads, candidate_generator_);
    
    const auto candidates = unique_left_align(candidate_generator_.generate_candidates(candidate_region), reference_);
    
    candidate_generator_.clear();
    
    //debug::print_candidates(candidates);
    
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
            break;
        }
        
        std::cout << "active region is " << active_region << '\n';
        //std::cout << "haplotype region is " << haplotypes.front().get_region() << '\n';
        
        const auto active_reads = copy_overlapped(reads, active_region);
        
        std::cout << "there are " << haplotypes.size() << " haplotypes" << '\n';
        //std::cout << "there are " << count_reads(active_reads) << " reads" << '\n';
        
        resume_timer(likelihood_timer);
        haplotype_likelihoods.populate(active_reads, haplotypes,
                                        flank_state(haplotypes, active_region, candidates));
        pause_timer(likelihood_timer);
        
        resume_timer(haplotype_fitler_timer);
        auto removed_haplotypes = filter_n_haplotypes(haplotypes, samples,
                                                      haplotype_likelihoods, max_haplotypes_);
        pause_timer(haplotype_fitler_timer);
        
        if (haplotypes.empty()) {
            //debug::print_read_haplotype_liklihoods(removed_haplotypes, active_reads, haplotype_likelihoods);
            // This can only happen if all haplotypes have equal likelihood.
            // TODO: is there anything else we can do?
            generator.remove_haplotypes(removed_haplotypes);
            continue;
        }
        
        resume_timer(likelihood_timer);
        haplotype_likelihoods.erase(removed_haplotypes);
        pause_timer(likelihood_timer);
        
        // Compute haplotype priors after likelihood filtering as prior model may use
        // interdependencies between haplotypes.
        resume_timer(prior_model_timer);
        auto haplotype_priors = haplotype_prior_model_->compute_maximum_entropy_haplotype_set(haplotypes);
        pause_timer(prior_model_timer);
        
        resume_timer(haplotype_generation_timer);
        generator.keep_haplotypes(haplotypes);
        generator.remove_haplotypes(removed_haplotypes);
        removed_haplotypes.clear();
        pause_timer(haplotype_generation_timer);
        
        //std::cout << "there are " << haplotypes.size() << " final haplotypes" << '\n';
        debug::print_read_haplotype_liklihoods(haplotypes, active_reads, haplotype_likelihoods, -1);
        
        resume_timer(latent_timer);
        const auto caller_latents = infer_latents(samples, haplotypes, haplotype_priors,
                                                  haplotype_likelihoods);
        pause_timer(latent_timer);
        
        haplotype_priors.clear();
        haplotype_likelihoods.clear();
        
        //debug::print_haplotype_posteriors(*caller_latents->get_haplotype_posteriors());
        
        resume_timer(phasing_timer);
        const auto phase_set = phaser.try_phase(haplotypes, *caller_latents->get_genotype_posteriors(),
                                                candidates);
        pause_timer(phasing_timer);
        
        auto unphased_active_region = active_region;
        
        if (overlaps(active_region, call_region) && phase_set) {
            assert(!is_empty_region(phase_set->region));
            
            //std::cout << "phased region is " << phase_set->region << '\n';
            
            auto active_candidates = copy_overlapped(candidates, phase_set->region);
            
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
            
            //std::cout << "removing " << removable_haplotypes.size() << " haplotypes" << '\n';
            
            resume_timer(haplotype_generation_timer);
            generator.remove_haplotypes(removable_haplotypes);
            next_active_region = generator.tell_next_active_region();
            pause_timer(haplotype_generation_timer);
        }
        
        if (overlaps(active_region, call_region) && has_moved_forward(next_active_region, active_region)) {
            auto passed_region = left_overhang_region(active_region, next_active_region);
            
            auto uncalled_region = overlapped_region(active_region, passed_region);
            
            if (phase_set && ends_before(phase_set->region, passed_region)) {
                uncalled_region = right_overhang_region(passed_region, phase_set->region);
            }
            
            const auto active_candidates = copy_overlapped(candidates, uncalled_region);
            
            resume_timer(phasing_timer);
            //std::cout << "force phasing " << uncalled_region << '\n';
            const auto phasings = phaser.force_phase(haplotypes,
                                                     *caller_latents->get_genotype_posteriors(),
                                                     active_candidates);
            pause_timer(phasing_timer);
            
            haplotypes.clear();
            haplotypes.shrink_to_fit();
            
            resume_timer(allele_generator_timer);
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
        
        debug::print_progress(active_region, call_region, 10000);
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
        if (p.second < 1e-06) {
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
    void print_candidates(const std::vector<Variant>& candidates, bool number_only)
    {
        if (candidates.empty()) {
            std::cout << "there are no candidates" << '\n';
        } else {
            std::cout << "found " << candidates.size() << " candidates: " << '\n';
            if (!number_only) {
                for (const auto& c : candidates) std::cout << c << '\n';
            }
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
    
    double percent_completed(const GenomicRegion& active_region,
                             const GenomicRegion& call_region)
    {
        const auto num_bases_processed = active_region.get_end() - call_region.get_begin();
        const auto total_bases = region_size(call_region);
        return 100 * static_cast<double>(num_bases_processed) / total_bases;
    }
    
    void print_progress(const GenomicRegion& last_active_region,
                        const GenomicRegion& call_region,
                        GenomicRegion::SizeType step_size)
    {
        auto& log = Logging::logger::get();
        
        BOOST_LOG_SEV(log, Logging::logging::trivial::info)
            << last_active_region.get_contig_name() << ':' << last_active_region.get_begin()
            << "\t"
            << percent_completed(last_active_region, call_region) << "%";
    }
    
    template <typename Map>
    void print_haplotype_posteriors(const Map& haplotype_posteriors, std::size_t n)
    {
        auto m = std::min(haplotype_posteriors.size(), n);
        
        std::cout << "printing top " << m << " haplotype posteriors" << std::endl;
        
        std::vector<std::pair<Haplotype, double>> v {};
        v.reserve(haplotype_posteriors.size());
        
        std::copy(std::cbegin(haplotype_posteriors), std::cend(haplotype_posteriors),
                  std::back_inserter(v));
        
        std::sort(std::begin(v), std::end(v),
                  [] (const auto& lhs, const auto& rhs) {
                      return lhs.second > rhs.second;
                  });
        
        for (unsigned i {}; i < m; ++i) {
            print_variant_alleles(v[i].first);
            std::cout << " " << std::setprecision(10) << v[i].second << std::endl;
        }
    }
}
} // namespace Octopus
