// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "caller.hpp"

#include <algorithm>
#include <utility>
#include <tuple>
#include <iterator>
#include <stdexcept>
#include <cassert>
#include <iostream>

#include "concepts/mappable.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/read_stats.hpp"
#include "utils/maths.hpp"
#include "core/models/haplotype_likelihood_model.hpp"
#include "core/tools/haplotype_filter.hpp"
#include "utils/call.hpp"
#include "utils/variant_call.hpp"
#include "utils/reference_call.hpp"

#include "timers.hpp"

namespace octopus {

// public methods

Caller::Caller(Components&& components, Parameters parameters)
: reference_ {components.reference}
, samples_ {components.read_pipe.get().samples()}
, debug_log_ {}
, trace_log_ {}
, read_pipe_ {components.read_pipe}
, candidate_generator_ {std::move(components.candidate_generator)}
, haplotype_generator_builder_ {std::move(components.haplotype_generator_builder)}
, phaser_ {std::move(components.phaser)}
, parameters_ {std::move(parameters)}
{
    if (parameters_.max_haplotypes == 0) {
        throw std::logic_error {"Caller: max haplotypes must be > 0"};
    }
    
    if (DEBUG_MODE) {
        debug_log_ = logging::DebugLogger {};
    }
    if (TRACE_MODE) {
        trace_log_ = logging::TraceLogger {};
    }
}

std::string Caller::name() const
{
    return do_name();
}

Caller::CallTypeSet Caller::call_types() const
{
    return do_call_types();
}

namespace debug {
    template <typename S>
    void print_left_aligned_candidates(S&& stream, const std::vector<Variant>& raw_candidates,
                                       const ReferenceGenome& reference);
    void print_left_aligned_candidates(const std::vector<Variant>& raw_candidates,
                                       const ReferenceGenome& reference);
    
    template <typename S>
    void print_final_candidates(S&& stream, const MappableFlatSet<Variant>& candidates,
                                bool number_only = false);
    void print_final_candidates(const MappableFlatSet<Variant>& candidates,
                                bool number_only = false);
    
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
    
    enum class Resolution {
        sequence, alleles, variantAlleles,
        sequenceAndAlleles, SequenceAndVariantAlleles
    };
    
    template <typename S>
    void print_haplotypes(S&& stream, const std::vector<Haplotype>& haplotypes,
                          Resolution resolution = Resolution::sequenceAndAlleles);
    void print_haplotypes(const std::vector<Haplotype>& haplotypes,
                          Resolution resolution = Resolution::sequenceAndAlleles);
    
    template <typename S, typename Map>
    void print_haplotype_posteriors(S&& stream, const Map& haplotype_posteriors, std::size_t n = 5);
    template <typename Map>
    void print_haplotype_posteriors(const Map& haplotype_posteriors, std::size_t n = 5);
    
    auto find_read(const std::string& region, const std::string& cigar_str,
                   const ReadContainer& reads);
    
    auto find_read(const SampleName& sample, const std::string& region,
                   const std::string& cigar_str, const ReadMap& reads);
    
    auto find_first_read(const std::string& region, const std::string& cigar_str,
                         const ReadMap& reads);
    
    double calculate_likelihood(const Haplotype& haplotype, const AlignedRead& read,
                                const HaplotypeLikelihoodModel::FlankState flank_state);
    
    void run_likelihood_calculation(const std::string& haplotype_str,
                                    const std::string& haplotype_region,
                                    const std::string& active_region,
                                    const std::string& read_region,
                                    const std::string& cigar_str,
                                    const ReadMap& reads,
                                    const MappableFlatSet<Variant>& candidates,
                                    const ReferenceGenome& reference);
}

bool is_boundry_insertion(const Variant& variant, const GenomicRegion& region)
{
    return is_empty_region(variant) && ends_equal(variant, region);
}

auto copy_contained_to_vector(const MappableFlatSet<Variant>& candidates,
                              const GenomicRegion& region,
                              const bool remove_rhs_boundry_insertions = false)
{
    auto contained = contained_range(candidates, region);
    if (remove_rhs_boundry_insertions) {
        while (!empty(contained) && is_boundry_insertion(contained.back(), region)) {
            contained.advance_end(-1);
        }
    }
    return std::vector<Variant> {std::cbegin(contained), std::cend(contained)};
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
        candidates.erase_overlapped(passed_region);
    }
}

template <typename Container>
void remove_duplicates(Container& haplotypes, const ReferenceGenome& reference)
{
    const auto n = unique_least_complex(haplotypes, Haplotype {haplotype_region(haplotypes), reference});
    if (DEBUG_MODE) {
        logging::DebugLogger log {};
        stream(log) << n << " duplicate haplotypes were removed";
    }
}

bool all_empty(const ReadMap& reads)
{
    return std::all_of(std::cbegin(reads), std::cend(reads),
                       [] (const auto& p) { return p.second.empty(); });
}

auto calculate_candidate_region(const GenomicRegion& call_region, const ReadMap& reads,
                                const VariantGenerator& candidate_generator)
{
    if (!candidate_generator.requires_reads()) return call_region;
    return all_empty(reads) ? call_region : encompassing_region(reads);
}

bool has_passed(const GenomicRegion& next_active_region, const GenomicRegion& active_region)
{
    return is_after(next_active_region, active_region) && active_region != next_active_region;
}

// Wrap the pointer so can use mappable algorithms
struct CallWrapper : public Mappable<CallWrapper>
{
    CallWrapper(std::unique_ptr<Call> call) : call {std::move(call) } {}
    operator const std::unique_ptr<Call>&() const noexcept { return call; }
    operator std::unique_ptr<Call>&() noexcept { return call; }
    std::unique_ptr<Call>::pointer operator->() const noexcept { return call.get(); };
    std::unique_ptr<Call> call;
    const GenomicRegion& mapped_region() const noexcept { return call->mapped_region(); }
};

template <typename T>
auto wrap(std::vector<std::unique_ptr<T>>&& calls)
{
    return std::vector<CallWrapper> {
        std::make_move_iterator(std::begin(calls)), std::make_move_iterator(std::end(calls))
    };
}

auto unwrap(std::vector<CallWrapper>&& calls)
{
    std::vector<std::unique_ptr<Call>> result {};
    result.reserve(calls.size());
    
    std::transform(std::make_move_iterator(std::begin(calls)),
                   std::make_move_iterator(std::end(calls)),
                   std::back_inserter(result),
                   [] (CallWrapper&& wrapped_call) -> std::unique_ptr<Call> {
                       return std::move(wrapped_call.call);
                   });
    
    calls.clear();
    calls.shrink_to_fit();
    
    return result;
}

void set_phase(const SampleName& sample, const Phaser::PhaseSet::PhaseRegion& phase,
               const std::vector<GenomicRegion>& call_regions, CallWrapper& call)
{
    const auto overlapped = overlap_range(call_regions, phase.mapped_region(),
                                          BidirectionallySortedTag {});
    
    if (!overlapped.empty()) {
        call->set_phase(sample, Call::PhaseCall {encompassing_region(overlapped.front(), call), phase.score});
    } else {
        call->set_phase(sample, Call::PhaseCall {mapped_region(call), phase.score});
    }
}

void set_phasing(std::vector<CallWrapper>& calls, const Phaser::PhaseSet& phase_set,
                 const GenomicRegion& calling_region)
{
    if (!calls.empty()) {
        const auto call_regions = extract_regions(calls);
        for (auto& call : calls) {
            const auto& call_region = mapped_region(call);
            
            for (const auto& p : phase_set.phase_regions) {
                const auto phase = find_phase_region(p.second, call_region);
                
                if (phase && overlaps(calling_region, phase->get().region)) {
                    if (begins_before(phase->get().region, calling_region)) {
                        const auto output_call_region = *overlapped_region(calling_region, phase->get().region);
                        const auto output_calls = overlap_range(call_regions, output_call_region);
                        
                        if (!output_calls.empty()) {
                            const Phaser::PhaseSet::PhaseRegion clipped_phase {
                                expand_lhs(phase->get().region, begin_distance(output_calls.front(), phase->get().region)),
                                phase->get().score
                            };
                            set_phase(p.first, clipped_phase, call_regions, call);
                        }
                    } else {
                        set_phase(p.first, *phase, call_regions, call);
                    }
                }
            }
        }
    }
}

namespace {
    auto mapped_region(const VcfRecord& record)
    {
        using SizeType = GenomicRegion::Size;
        const auto begin = record.pos() - 1;
        return GenomicRegion {record.chrom(), begin, begin + static_cast<SizeType>(record.ref().size())};
    }
} // namespace

void erase_calls_outside_region(std::vector<VcfRecord>& calls, const GenomicRegion& region)
{
    auto it = std::remove_if(std::begin(calls), std::end(calls),
                             [&region] (const auto& call) {
                                 return !overlaps(mapped_region(call), region);
                             });
    calls.erase(it, std::end(calls));
}

void merge(std::vector<CallWrapper>&& src, std::deque<VcfRecord>& dst,
           const VcfRecordFactory& factory, const GenomicRegion& call_region)
{
    using std::begin; using std::end;
    if (src.empty()) return;
    auto new_records = factory.make(unwrap(std::move(src)));
    erase_calls_outside_region(new_records, call_region);
    const auto it = dst.insert(end(dst),
                               std::make_move_iterator(begin(new_records)),
                               std::make_move_iterator(end(new_records)));
    std::inplace_merge(begin(dst), it, end(dst),
                       [] (const auto& lhs, const auto& rhs) {
                           return mapped_region(lhs) < mapped_region(rhs);
                       });
    // sometimes duplicates are called on active region boundries
    const auto it2 = std::unique(begin(dst), end(dst),
                                 [] (const auto& lhs, const auto& rhs) {
                                     return lhs.pos() == rhs.pos() && lhs.ref() == rhs.ref()
                                            && lhs.alt() == rhs.alt();
                                 });
    dst.erase(it2, end(dst));
}

std::deque<VcfRecord> Caller::call(const GenomicRegion& call_region, ProgressMeter& progress_meter) const
{
    // TODO: Needs refactoring into smaller methods!
    resume(init_timer);
    
    ReadMap reads;
    std::deque<VcfRecord> result {};
    
    if (candidate_generator_.requires_reads()) {
        reads = read_pipe_.get().fetch_reads(call_region);
        add_reads(reads, candidate_generator_);
        if (!refcalls_requested() && all_empty(reads)) {
            if (debug_log_) stream(*debug_log_) << "Stopping early as no reads found in call region";
            return result;
        }
        if (debug_log_) stream(*debug_log_) << "There are " << count_reads(reads) << " reads";
    }
    
    const auto candidate_region = calculate_candidate_region(call_region, reads, candidate_generator_);
    auto candidates = generate_candidate_variants(candidate_region);
    
    if (debug_log_) debug::print_final_candidates(stream(*debug_log_), candidates);
    
    if (!refcalls_requested() && candidates.empty()) {
        progress_meter.log_completed(call_region);
        return result;
    }
    if (!candidate_generator_.requires_reads()) {
        // as we didn't fetch them earlier
        reads = read_pipe_.get().fetch_reads(extract_regions(candidates));
    }
    
    auto haplotype_generator   = make_haplotype_generator(candidates, reads);
    auto haplotype_likelihoods = make_haplotype_likelihood_cache();
    const auto record_factory  = make_record_factory(reads);
    std::vector<Haplotype> haplotypes {}, next_haplotypes {};
    GenomicRegion active_region;
    boost::optional<GenomicRegion> next_active_region {};
    auto completed_region = head_region(call_region);
    
    pause(init_timer);
    
    while (true) {
        if (next_active_region) {
            haplotypes = std::move(next_haplotypes);
            active_region = std::move(*next_active_region);
            next_active_region = boost::none;
        } else {
            try {
                resume(haplotype_generation_timer);
                std::tie(haplotypes, active_region, std::ignore) = haplotype_generator.generate();
                pause(haplotype_generation_timer);
            } catch (const HaplotypeGenerator::HaplotypeOverflow& e) {
                logging::WarningLogger wlog {};
                stream(wlog) << "Skipping region " << e.region() << " as there are too many haplotypes";
                haplotype_generator.clear_progress();
                haplotype_likelihoods.clear();
                progress_meter.log_completed(e.region());
                continue;
            }
        }
        
        if (debug_log_) stream(*debug_log_) << "Active region is " << active_region;
        
        if (is_after(active_region, call_region) || haplotypes.empty()) {
            if (debug_log_) {
                if (haplotypes.empty()) {
                    stream(*debug_log_) << "No haplotypes were generated in the active region";
                } else {
                    stream(*debug_log_) << "Generated " << haplotypes.size()
                                        << " haplotypes but active region is after call region";
                }
            }
            progress_meter.log_completed(active_region);
            break;
        }
        
        const auto active_reads = copy_overlapped(reads, active_region);
        
        if (debug_log_) stream(*debug_log_) << "There are " << haplotypes.size() << " initial haplotypes";
        
        if (!refcalls_requested() && !has_coverage(active_reads)) {
            if (debug_log_) stream(*debug_log_) << "Skipping active region as there are no active reads";
            continue;
        }
        if (debug_log_) {
            stream(*debug_log_) << "There are " << count_reads(active_reads) << " active reads";
        }
        
        remove_duplicates(haplotypes, reference_);
        
        try {
            resume(haplotype_likelihood_timer);
            populate(haplotype_likelihoods, active_region, haplotypes, candidates, active_reads);
            pause(haplotype_likelihood_timer);
        } catch(const HaplotypeLikelihoodModel::ShortHaplotypeError& e) {
            if (debug_log_) {
                stream(*debug_log_) << "Skipping " << active_region << " as a haplotype was too short by "
                                    << e.required_extension() << "bp";
            }
            // TODO: we could force HaplotypeGenerator to extend the current set of haplotypes
            // and retry
            haplotype_generator.clear_progress();
            haplotype_likelihoods.clear();
            continue;
        }
        
        auto removed_haplotypes = filter(haplotypes, haplotype_likelihoods);
        
        if (haplotypes.empty()) {
            // This can only happen if all haplotypes have equal likelihood
            haplotype_generator.clear_progress();
            haplotype_likelihoods.clear();
            continue;
        }
        if (haplotypes.capacity() > 2 * haplotypes.size()) {
            haplotypes.shrink_to_fit();
        }
        
        resume(haplotype_likelihood_timer);
        haplotype_likelihoods.erase(removed_haplotypes);
        pause(haplotype_likelihood_timer);
        
        resume(haplotype_generation_timer);
        auto has_removal_impact = haplotype_generator.removal_has_impact();
        if (has_removal_impact) {
            try {
                haplotype_generator.remove(removed_haplotypes);
                haplotype_generator.collapse(haplotypes);
            } catch (const std::runtime_error& e) {
                std::cout << "Exception: " << e.what() << ". Encountered in active region "
                          << active_region << " whilst calling " << call_region << std::endl;
                throw;
            }
        } else {
            haplotype_generator.clear_progress();
        }
        pause(haplotype_generation_timer);
        removed_haplotypes.clear();
        removed_haplotypes.shrink_to_fit();
        
        if (debug_log_) stream(*debug_log_) << "There are " << haplotypes.size() << " final haplotypes";
        
        resume(latent_timer);
        const auto caller_latents = this->infer_latents(haplotypes, haplotype_likelihoods);
        pause(latent_timer);
        
        if (trace_log_) {
            debug::print_haplotype_posteriors(stream(*trace_log_), *caller_latents->haplotype_posteriors(), -1);
        } else if (debug_log_) {
            debug::print_haplotype_posteriors(stream(*debug_log_), *caller_latents->haplotype_posteriors());
        }
        
        resume(phasing_timer);
        const auto phase_set = phaser_.try_phase(haplotypes, *caller_latents->genotype_posteriors(),
                                                 copy_contained_to_vector(candidates, haplotype_region(haplotypes)));
        pause(phasing_timer);
        
        if (debug_log_) {
            if (phase_set) {
                debug::print_phase_sets(stream(*debug_log_), *phase_set);
            } else {
                *debug_log_ << "No partial phasings found";
            }
        }
        
        if (phase_set && overlaps(active_region, call_region)) {
            assert(!is_empty(phase_set->region));
            
            if (debug_log_) stream(*debug_log_) << "Phased region is " << phase_set->region;
            
            const auto active_candidates = copy_contained_to_vector(candidates, phase_set->region);
            
            if (!active_candidates.empty()) {
                auto variant_calls = wrap(call_variants(active_candidates, *caller_latents));
                if (!variant_calls.empty()) {
                    if (parameters_.allow_model_filtering) {
                        const auto mp = calculate_model_posterior(haplotypes, haplotype_likelihoods,
                                                                  *caller_latents);
                        if (mp) {
                            for (auto& call : variant_calls) {
                                call->set_model_posterior(probability_to_phred(1 - *mp));
                            }
                        }
                    }
                    set_phasing(variant_calls, *phase_set, call_region);
                    merge(std::move(variant_calls), result, record_factory, call_region);
                }
            }
            
            active_region = right_overhang_region(active_region, phase_set->region);
            haplotype_generator.jump(active_region);
        } else {
            resume(haplotype_generation_timer);
            if (has_removal_impact) { // if there was no impact before then there can't be now either
                has_removal_impact = haplotype_generator.removal_has_impact();
            }
            if (has_removal_impact) {
                const auto max_to_remove = haplotype_generator.max_removal_impact();
                auto removable_haplotypes = get_removable_haplotypes(haplotypes,
                                                                     haplotype_likelihoods,
                                                                     *caller_latents->haplotype_posteriors(),
                                                                     max_to_remove);
                if (debug_log_) {
                    stream(*debug_log_) << "Discarding " << removable_haplotypes.size()
                                    << " haplotypes with low posterior support";
                }
                haplotype_generator.remove(removable_haplotypes);
            }
            pause(haplotype_generation_timer);
        }
        
        if (!parameters_.allow_model_filtering) {
            haplotype_likelihoods.clear();
        }
        
        resume(haplotype_generation_timer);
        bool last_pass;
        try {
            resume(haplotype_generation_timer);
            std::tie(next_haplotypes, next_active_region, last_pass) = haplotype_generator.generate();
            pause(haplotype_generation_timer);
        } catch (const HaplotypeGenerator::HaplotypeOverflow& e) {
            logging::WarningLogger wlog {};
            stream(wlog) << "Skipping region " << e.region() << " as there are too many haplotypes";
            haplotype_generator.clear_progress();
            next_active_region = tail_region(e.region());
        }
        pause(haplotype_generation_timer);
        
        if (last_pass && begins_before(active_region, *next_active_region) && overlaps(active_region, call_region)) {
            auto passed_region   = left_overhang_region(active_region, *next_active_region);
            auto uncalled_region = *overlapped_region(active_region, passed_region);
            
            if (phase_set && ends_before(phase_set->region, passed_region)) {
                uncalled_region = right_overhang_region(passed_region, phase_set->region);
            }
            
            auto active_candidates = copy_contained_to_vector(candidates, uncalled_region,
                                                              are_adjacent(uncalled_region, *next_active_region));
            
            std::vector<GenomicRegion> called_regions;
            
            if (!active_candidates.empty()) {
                if (debug_log_) stream(*debug_log_) << "Calling variants in region " << uncalled_region;
    
                resume(calling_timer);
                auto variant_calls = wrap(call_variants(active_candidates, *caller_latents));
                pause(calling_timer);
                
                if (!variant_calls.empty()) {
                    if (parameters_.allow_model_filtering) {
                        const auto mp = calculate_model_posterior(haplotypes, haplotype_likelihoods, *caller_latents);
                        if (mp) {
                            for (auto& call : variant_calls) {
                                call->set_model_posterior(probability_to_phred(1 - *mp));
                            }
                        }
                    }
                    called_regions = extract_covered_regions(variant_calls);
    
                    resume(phasing_timer);
                    const auto phase = phaser_.force_phase(haplotypes,
                                                           *caller_latents->genotype_posteriors(),
                                                           active_candidates);
                    pause(phasing_timer);
                    
                    if (debug_log_) debug::print_phase_sets(stream(*debug_log_), phase);
    
                    resume(misc_timer[3]);
                    set_phasing(variant_calls, phase, call_region);
                    merge(std::move(variant_calls), result, record_factory, call_region);
                    pause(misc_timer[3]);
                }
            }
            
            if (refcalls_requested()) {
                auto alleles = generate_candidate_reference_alleles(uncalled_region, active_candidates, called_regions);
                auto reference_calls = wrap(this->call_reference(alleles, *caller_latents, reads));
                merge(std::move(reference_calls), result, record_factory, call_region);
            }
            
            completed_region = encompassing_region(completed_region, passed_region);
        }
        
        haplotype_likelihoods.clear();
        progress_meter.log_completed(completed_region);
    }
    
    return result;
}
    
std::vector<VcfRecord> Caller::regenotype(const std::vector<Variant>& variants, ProgressMeter& progress_meter) const
{
    return {}; // TODO
}

// private methods

bool Caller::refcalls_requested() const noexcept
{
    return parameters_.refcall_type != RefCallType::none;
}

MappableFlatSet<Variant> Caller::generate_candidate_variants(const GenomicRegion& region) const
{
    if (debug_log_) stream(*debug_log_) << "Generating candidate variants in region " << region;
    auto raw_candidates = candidate_generator_.generate(region);
    if (debug_log_) debug::print_left_aligned_candidates(stream(*debug_log_), raw_candidates, reference_);
    auto final_candidates = unique_left_align(std::move(raw_candidates), reference_);
    candidate_generator_.clear();
    return MappableFlatSet<Variant> {
        std::make_move_iterator(std::begin(final_candidates)),
        std::make_move_iterator(std::end(final_candidates))
    };
}

HaplotypeGenerator Caller::make_haplotype_generator(const MappableFlatSet<Variant>& candidates,
                                                    const ReadMap& reads) const
{
    return haplotype_generator_builder_.build(reference_, candidates, reads);
}

HaplotypeLikelihoodCache Caller::make_haplotype_likelihood_cache() const
{
    return HaplotypeLikelihoodCache {parameters_.max_haplotypes, samples_};
}

VcfRecordFactory Caller::make_record_factory(const ReadMap& reads) const
{
    return VcfRecordFactory {reference_, reads, samples_, parameters_.call_sites_only};
}

auto calculate_flank_regions(const GenomicRegion& haplotype_region,
                             const GenomicRegion& active_region,
                             const MappableFlatSet<Variant>& candidates)
{
    auto lhs_flank = left_overhang_region(haplotype_region, active_region);
    auto rhs_flank = right_overhang_region(haplotype_region, active_region);
    
    const auto active_candidates = contained_range(candidates, active_region);
    
    assert(!active_candidates.empty());
    
    if (is_empty_region(*leftmost_mappable(active_candidates)) && !is_empty(lhs_flank)) {
        lhs_flank = expand_rhs(lhs_flank, -1); // stops boundry insertions being inactive
    }
    
    const auto lhs_inactive_candidates = contained_range(candidates, lhs_flank);
    if (lhs_inactive_candidates.empty()) {
        lhs_flank = head_region(lhs_flank);
    } else {
        lhs_flank = closed_region(lhs_flank, rightmost_region(lhs_inactive_candidates));
    }
    
    if (is_empty_region(*rightmost_mappable(active_candidates)) && !is_empty(rhs_flank)) {
        rhs_flank = expand_lhs(rhs_flank, -1); // stops boundry insertions being inactive
    }
    
    const auto rhs_inactive_candidates = contained_range(candidates, rhs_flank);
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
    const auto flank_regions = calculate_flank_regions(haplotype_region(haplotypes), active_region,
                                                       candidates);
    return HaplotypeLikelihoodCache::FlankState {
        size(flank_regions.first), size(flank_regions.second)
    };
}

void Caller::populate(HaplotypeLikelihoodCache& haplotype_likelihoods,
                      const GenomicRegion& active_region,
                      const std::vector<Haplotype>& haplotypes,
                      const MappableFlatSet<Variant>& candidates,
                      const ReadMap& active_reads) const
{
    assert(haplotype_likelihoods.is_empty());
    
    boost::optional<HaplotypeLikelihoodCache::FlankState> flank_state {};
    
    if (debug_log_) {
        stream(*debug_log_) << "Calculating likelihoods for " << haplotypes.size() << " haplotypes";
        debug::print_active_candidates(stream(*debug_log_), candidates, active_region);
        stream(*debug_log_) << "Haplotype region is " << haplotype_region(haplotypes);
    }
    
    if (parameters_.allow_inactive_flank_scoring) {
        flank_state = calculate_flank_state(haplotypes, active_region, candidates);
        if (debug_log_) {
            debug::print_inactive_flanking_candidates(stream(*debug_log_), candidates, active_region,
                                                      haplotype_region(haplotypes));
        }
    }
    
    haplotype_likelihoods.populate(active_reads, haplotypes, std::move(flank_state));
    
    if (trace_log_) {
        debug::print_read_haplotype_likelihoods(stream(*trace_log_), haplotypes, active_reads,
                                               haplotype_likelihoods, -1);
    }
}

std::vector<Haplotype>
Caller::filter(std::vector<Haplotype>& haplotypes, const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    auto removed_haplotypes = filter_to_n(haplotypes, samples_, haplotype_likelihoods,
                                          parameters_.max_haplotypes);
    
    if (debug_log_) {
        if (haplotypes.empty()) {
            *debug_log_ << "Filtered all haplotypes";
        } else {
            stream(*debug_log_) << "Filtered " << removed_haplotypes.size() << " haplotypes";
        }
    }
    if (trace_log_) {
        stream(*trace_log_) << "Filtered " << removed_haplotypes.size() << " haplotypes:";
        debug::print_haplotypes(stream(*trace_log_), removed_haplotypes,
                                debug::Resolution::variantAlleles);
    }
    
    return removed_haplotypes;
}

std::vector<std::reference_wrapper<const Haplotype>>
Caller::get_removable_haplotypes(const std::vector<Haplotype>& haplotypes,
                                 const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                 const Caller::Latents::HaplotypeProbabilityMap& haplotype_posteriors,
                                 const unsigned max_to_remove) const
{
    return extract_removable(haplotypes, haplotype_posteriors, samples_, haplotype_likelihoods,
                             max_to_remove, parameters_.haplotype_extension_threshold.probability_false());
}

bool Caller::done_calling(const GenomicRegion& region) const noexcept
{
    return is_empty(region);
}

std::vector<Allele>
Caller::generate_callable_alleles(const GenomicRegion& region, const std::vector<Variant>& candidates) const
{
    using std::begin; using std::end; using std::make_move_iterator; using std::back_inserter;
    
    auto overlapped_candidates = copy_overlapped(candidates, region);
    
    if (is_empty(region) && overlapped_candidates.empty()) return {};
    
    if (overlapped_candidates.empty()) {
        switch (parameters_.refcall_type) {
            case RefCallType::positional:
                return make_positional_reference_alleles(region, reference_);
            case RefCallType::blocked:
                return std::vector<Allele> {make_reference_allele(region, reference_)};
            default:
                return {};
        }
    }
    
    auto variant_alleles = decompose(overlapped_candidates);
    
    if (parameters_.refcall_type == RefCallType::none) return variant_alleles;
    
    auto covered_regions   = extract_covered_regions(overlapped_candidates);
    auto uncovered_regions = extract_intervening_regions(covered_regions, region);
    
    std::vector<Allele> result {};
    
    if (parameters_.refcall_type == Caller::RefCallType::blocked) {
        auto reference_alleles = make_reference_alleles(uncovered_regions, reference_);
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
                auto alleles = make_positional_reference_alleles(*uncovered_itr, reference_);
                result.insert(end(result),
                              make_move_iterator(begin(alleles)),
                              make_move_iterator(end(alleles)));
                std::advance(uncovered_itr, 1);
            }
            result.push_back(std::move(variant_allele));
        }
        
        if (uncovered_itr != uncovered_end) {
            auto alleles = make_positional_reference_alleles(*uncovered_itr, reference_);
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
                   const Caller::RefCallType refcall_type)
{
    if (refcall_type == Caller::RefCallType::blocked && !alleles.empty()
        && are_adjacent(alleles.back(), allele)) {
        alleles.back() = Allele {encompassing_region(alleles.back(), allele),
            alleles.back().sequence() + allele.sequence()};
    } else {
        alleles.push_back(allele);
    }
}

// TODO: we should catch the case where an insertion has been called and push the refcall
// block up a position, otherwise the returned reference allele (block) will never be called.
std::vector<Allele>
Caller::generate_candidate_reference_alleles(const GenomicRegion& region,
                                             const std::vector<Variant>& candidates,
                                             const std::vector<GenomicRegion>& called_regions) const
{
    using std::cbegin; using std::cend;
    
    auto callable_alleles = generate_callable_alleles(region, candidates);
    
    if (callable_alleles.empty() || parameters_.refcall_type == RefCallType::none) return {};
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
            append_allele(result, *allele_itr, parameters_.refcall_type);
            std::copy(std::next(allele_itr), allele_end_itr, std::back_inserter(result));
            break;
        }
        if (called_itr == called_end_itr) {
            append_allele(result, *allele_itr, parameters_.refcall_type);
            if (begins_before(*allele_itr, *candidate_itr)) {
                ++allele_itr;
            } else {
                allele_itr = find_next(allele_itr, allele_end_itr, *candidate_itr);
                ++candidate_itr;
            }
        } else {
            if (is_same_region(*called_itr, *candidate_itr)) { // called candidate
                while (is_before(*allele_itr, *called_itr)) {
                    append_allele(result, *allele_itr, parameters_.refcall_type);
                    ++allele_itr;
                }
                allele_itr = find_next(allele_itr, allele_end_itr, *candidate_itr);
                ++candidate_itr;
                ++called_itr;
            } else if (begins_before(*called_itr, *candidate_itr)) { // parsimonised called candidate
                if (!overlaps(*allele_itr, *called_itr)) {
                    append_allele(result, *allele_itr, parameters_.refcall_type);
                    ++allele_itr;
                } else {
                    if (begins_before(*allele_itr, *called_itr)) { // when variant has been left padded
                        append_allele(result, splice(*allele_itr, left_overhang_region(*allele_itr, *called_itr)),
                                      parameters_.refcall_type);
                    }
                    // skip contained alleles and candidates as they include called variants
                    allele_itr    = cend(contained_range(allele_itr, allele_end_itr, *called_itr)).base();
                    candidate_itr = cend(contained_range(candidate_itr, candidate_end_itr, *called_itr)).base();
                    ++called_itr;
                }
            } else {
                append_allele(result, *allele_itr, parameters_.refcall_type);
                if (begins_before(*allele_itr, *candidate_itr)) {
                    ++allele_itr;
                } else {
                    allele_itr = find_next(allele_itr, allele_end_itr, *candidate_itr);
                    ++candidate_itr;
                }
            }
        }
    }
    return result;
}

namespace debug {
    template <typename S>
    void print_left_aligned_candidates(S&& stream, const std::vector<Variant>& raw_candidates,
                                       const ReferenceGenome& reference)
    {
        std::deque<std::pair<Variant, Variant>> left_aligned {};
        
        for (const auto& raw_candidate : raw_candidates) {
            auto left_aligned_candidate = left_align(raw_candidate, reference);
            
            if (left_aligned_candidate != raw_candidate) {
                left_aligned.emplace_back(raw_candidate, std::move(left_aligned_candidate));
            }
        }
        
        if (left_aligned.empty()) {
            stream << "No candidates were left aligned" << '\n';
        } else {
            if (left_aligned.size() == 1) {
                stream << "1 candidate was left aligned:" << '\n';
            } else {
                stream << left_aligned.size() << " candidates were left aligned:" << '\n';
            }
            for (const auto& p : left_aligned) {
                stream << p.first << " to " << p.second << '\n';
            }
        }
    }
    
    void print_left_aligned_candidates(const std::vector<Variant>& raw_candidates,
                                       const ReferenceGenome& reference)
    {
        print_left_aligned_candidates(std::cout, raw_candidates, reference);
    }
    
    template <typename S>
    void print_final_candidates(S&& stream, const MappableFlatSet<Variant>& candidates, bool number_only)
    {
        if (candidates.empty()) {
            stream << "There are no final candidates" << '\n';
        } else {
            if (candidates.size() == 1) {
                stream << "There is 1 final candidate:" << '\n';
            } else {
                stream << "There are " << candidates.size() << " final candidates:" << '\n';
            }
            
            if (!number_only) {
                for (const auto& c : candidates) stream << c << '\n';
            }
        }
    }
    
    void print_final_candidates(const MappableFlatSet<Variant>& candidates, bool number_only)
    {
        print_final_candidates(std::cout, candidates, number_only);
    }
    
    template <typename S>
    void print_active_candidates(S&& stream, const MappableFlatSet<Variant>& candidates,
                                 const GenomicRegion& active_region, bool number_only)
    {
        const auto active_candidates = contained_range(candidates, active_region);
        
        if (active_candidates.empty()) {
            stream << "There are no active candidates" << '\n';
        } else {
            const auto num_active = size(active_candidates);
            
            if (num_active == 1) {
                stream << "There is 1 active candidate:" << '\n';
            } else {
                stream << "There are " << num_active << " active candidates:" << '\n';
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
        
        stream << "Haplotype flank regions are " << flanks.first << " and " << flanks.second << '\n';
        
        const auto lhs_inactive_candidates = contained_range(candidates, flanks.first);
        const auto rhs_inactive_candidates = contained_range(candidates, flanks.second);
        
        const auto num_lhs_inactives = size(lhs_inactive_candidates);
        const auto num_rhs_inactives = size(rhs_inactive_candidates);
        
        if (lhs_inactive_candidates.empty()) {
            if (rhs_inactive_candidates.empty()) {
                stream << "There are no inactive flanking candidates" << '\n';
                return;
            }
            
            stream << "There are no lhs inactive flanking candidates" << '\n';
            
            if (num_rhs_inactives == 1) {
                stream << "There is 1 possible rhs inactive flanking candidates: " << '\n';
            } else {
                stream << "There are " << num_rhs_inactives << " possible rhs inactive flanking candidates: " << '\n';
            }
            
            if (!number_only) {
                for (const auto& c : rhs_inactive_candidates) stream << c << '\n';
            }
            
            return;
        }
        
        if (num_lhs_inactives == 1) {
            stream << "There is 1 possible lhs inactive flanking candidates: " << '\n';
        } else {
            stream << "There are " << num_lhs_inactives << " possible lhs inactive flanking candidates: " << '\n';
        }
        
        if (!number_only) {
            for (const auto& c : lhs_inactive_candidates) stream << c << '\n';
        }
        
        if (rhs_inactive_candidates.empty()) {
            stream << "There are no rhs inactive flanking candidates" << '\n';
        } else {
            if (num_rhs_inactives == 1) {
                stream << "There is 1 possible rhs inactive flanking candidates: " << '\n';
            } else {
                stream << "There are " << num_rhs_inactives << " possible rhs inactive flanking candidates: " << '\n';
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
        print_inactive_flanking_candidates(std::cout, candidates, active_region, haplotype_region,
                                           number_only);
    }
    
    template <typename S>
    void print_haplotypes(S&& stream, const std::vector<Haplotype>& haplotypes,
                          Resolution resolution)
    {
        stream << "Printing " << haplotypes.size() << " haplotypes" << '\n';
        
        for (const auto& haplotype : haplotypes) {
            if (resolution == Resolution::sequence || resolution == Resolution::sequenceAndAlleles
                || resolution == Resolution::SequenceAndVariantAlleles) {
                stream << haplotype << '\n';
            }
            if (resolution == Resolution::alleles || resolution == Resolution::sequenceAndAlleles) {
                debug::print_alleles(stream, haplotype); stream << '\n';
            } else if (resolution != Resolution::sequence) {
                debug::print_variant_alleles(stream, haplotype);
                stream << '\n';
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
            stream << "Printing all haplotype posteriors" << '\n';
        } else {
            stream << "Printing top " << m << " haplotype posteriors" << '\n';
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
                          debug::print_variant_alleles(stream, p.first);
                          stream << " " << p.second << '\n';
                      });
    }
    
    template <typename Map>
    void print_haplotype_posteriors(const Map& haplotype_posteriors, std::size_t n)
    {
        print_haplotype_posteriors(std::cout, haplotype_posteriors, n);
    }
    
    auto find_read(const std::string& region, const std::string& cigar_str,
                   const ReadContainer& reads)
    {
        const auto cigar = parse_cigar(cigar_str);
        return std::find_if(std::cbegin(reads), std::cend(reads),
                            [&] (const AlignedRead& read) {
                                return read.cigar() == cigar
                                && to_string(mapped_region(read)) == region;
                            });
    }
    
    auto find_read(const SampleName& sample, const std::string& region,
                   const std::string& cigar_str, const ReadMap& reads)
    {
        return find_read(region, cigar_str, reads.at(sample));
    }
    
    auto find_first_read(const std::string& region, const std::string& cigar_str,
                         const ReadMap& reads)
    {
        for (const auto& p : reads) {
            const auto it = find_read(region, cigar_str, p.second);
            if (it != std::cend(p.second)) return it;
        }
        return std::cend(std::cbegin(reads)->second);
    }
    
    double calculate_likelihood(const Haplotype& haplotype, const AlignedRead& read,
                                HaplotypeLikelihoodModel::FlankState flank_state)
    {
        SampleName test_sample {"*test-sample*"};
        HaplotypeLikelihoodCache cache {1, {test_sample}};
        ReadContainer sample_reads {};
        sample_reads.emplace(read);
        ReadMap reads {};
        reads.emplace(test_sample, sample_reads);
        cache.populate(reads, {haplotype}, std::move(flank_state));
        return cache(test_sample, haplotype).front();
    }
    
    void run_likelihood_calculation(const std::string& haplotype_str,
                                    const std::string& haplotype_region_str,
                                    const std::string& active_region_str,
                                    const std::string& read_region_str,
                                    const std::string& cigar_str,
                                    const ReadMap& reads,
                                    const MappableFlatSet<Variant>& candidates,
                                    const ReferenceGenome& reference)
    {
//        auto haplotype = debug::make_haplotype(haplotype_str, haplotype_region_str, reference);
//        
//        std::cout << "Haplotype: " << haplotype << std::endl;
//        debug::print_variant_alleles(haplotype);
//        std::cout << std::endl;
//        
//        const auto active_region = parse_region(active_region_str, reference);
//        
//        auto flank_state = calculate_flank_state({haplotype}, active_region, candidates);
//        
//        std::cout << "Flank sizes: " << flank_state.lhs_flank << " " << flank_state.rhs_flank << std::endl;
//        
//        auto read = *find_first_read(read_region_str, cigar_str, reads);
//        
//        std::cout << "Read: " << read << std::endl;
//        
//        auto likelihood = calculate_likelihood(haplotype, read, flank_state);
//        
//        std::cout << "Likelihood = " << likelihood << std::endl;
    }
} // namespace debug
} // namespace octopus
