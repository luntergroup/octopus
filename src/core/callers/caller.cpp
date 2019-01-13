// Copyright (c) 2015-2018 Daniel Cooke
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
#include "core/types/calls/call.hpp"
#include "core/types/calls/call_wrapper.hpp"
#include "core/types/calls/variant_call.hpp"
#include "core/types/calls/reference_call.hpp"
#include "core/models/haplotype_likelihood_model.hpp"
#include "core/tools/haplotype_filter.hpp"
#include "core/tools/read_assigner.hpp"
#include "core/tools/read_realigner.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/read_stats.hpp"
#include "utils/maths.hpp"
#include "utils/append.hpp"

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
, likelihood_model_ {std::move(components.likelihood_model)}
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

unsigned Caller::min_callable_ploidy() const
{
    return do_min_callable_ploidy();
}

unsigned Caller::max_callable_ploidy() const
{
    return do_max_callable_ploidy();
}

namespace debug {

template <typename S>
void print_left_aligned_candidates(S&& stream, const std::vector<Variant>& raw_candidates,
                                   const ReferenceGenome& reference);
void print_left_aligned_candidates(const std::vector<Variant>& raw_candidates,
                                   const ReferenceGenome& reference);

template <typename S>
void print_final_candidates(S&& stream, const MappableFlatSet<Variant>& candidates, const GenomicRegion& region,
                            bool number_only = false);
void print_final_candidates(const MappableFlatSet<Variant>& candidates, const GenomicRegion& region,
                            bool number_only = false);

} // namespace debug

namespace {

bool all_empty(const ReadMap& reads)
{
    return std::all_of(std::cbegin(reads), std::cend(reads), [] (const auto& p) { return p.second.empty(); });
}

auto calculate_candidate_region(const GenomicRegion& call_region, const ReadMap& reads,
                                const ReferenceGenome& reference, const VariantGenerator& candidate_generator)
{
    if (!candidate_generator.requires_reads()) return call_region;
    const auto read_region = all_empty(reads) ? call_region : encompassing_region(reads);
    const auto contig_region = reference.contig_region(call_region.contig_name());
    if (right_overhangs(read_region, contig_region)) {
        return *overlapped_region(read_region, contig_region);
    } else {
        return read_region;
    }
}

auto mapped_region(const VcfRecord& record)
{
    using SizeType = GenomicRegion::Size;
    const auto begin = record.pos() - 1;
    return GenomicRegion {record.chrom(), begin, begin + static_cast<SizeType>(record.ref().size())};
}

void erase_calls_outside_region(std::vector<VcfRecord>& calls, const GenomicRegion& region)
{
    const auto overlapped = overlap_range(calls, region);
    calls.erase(overlapped.end().base(), std::cend(calls));
    calls.erase(std::cbegin(calls), overlapped.begin().base());
}

template <typename T>
auto to_vector(std::deque<T>&& values)
{
    using std::make_move_iterator;
    return std::vector<T> {make_move_iterator(std::begin(values)), make_move_iterator(std::end(values))};
}

auto convert_to_vcf(std::deque<CallWrapper>&& calls, const VcfRecordFactory& factory, const GenomicRegion& call_region)
{
    auto records = factory.make(to_vector(std::move(calls)));
    erase_calls_outside_region(records, call_region);
    std::deque<VcfRecord> result {};
    utils::append(std::move(records), result);
    return result;
}

} // namespace

std::deque<VcfRecord> Caller::call(const GenomicRegion& call_region, ProgressMeter& progress_meter) const
{
    ReadPipe::Report reads_report {};
    ReadMap reads;
    if (candidate_generator_.requires_reads()) {
        reads = read_pipe_.get().fetch_reads(expand(call_region, 100), reads_report);
        add_reads(reads, candidate_generator_);
        if (!refcalls_requested() && all_empty(reads)) {
            if (debug_log_) stream(*debug_log_) << "Stopping early as no reads found in call region " << call_region;
            return {};
        }
        if (debug_log_) stream(*debug_log_) << "Using " << count_reads(reads) << " reads in call region " << call_region;
    }
    const auto candidate_region = calculate_candidate_region(call_region, reads, reference_, candidate_generator_);
    auto candidates = generate_candidate_variants(candidate_region);
    if (debug_log_) debug::print_final_candidates(stream(*debug_log_), candidates, candidate_region);
    if (!refcalls_requested() && candidates.empty()) {
        progress_meter.log_completed(call_region);
        return {};
    }
    if (!candidate_generator_.requires_reads()) {
        // as we didn't fetch them earlier
        reads = read_pipe_.get().fetch_reads(call_region, reads_report);
    }
    auto calls = call_variants(call_region, candidates, reads, reads_report, progress_meter);
    candidates.clear();
    candidates.shrink_to_fit();
    progress_meter.log_completed(call_region);
    const auto record_factory = make_record_factory(reads);
    if (debug_log_) stream(*debug_log_) << "Converting " << calls.size() << " calls made in " << call_region << " to VCF";
    return convert_to_vcf(std::move(calls), record_factory, call_region);
}

std::vector<VcfRecord> Caller::regenotype(const std::vector<Variant>& variants, ProgressMeter& progress_meter) const
{
    return {}; // TODO
}

auto assign_and_realign(const std::vector<AlignedRead>& reads, const Genotype<Haplotype>& genotype)
{
    auto result = compute_haplotype_support(genotype, reads, {AssignmentConfig::AmbiguousAction::first});
    for (auto& p : result) realign_to_reference(p.second, p.first);
    return result;
}

// private methods

namespace debug {

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

} // namespace debug

namespace {

const auto& haplotype_region(const std::vector<Haplotype>& haplotypes)
{
    return mapped_region(haplotypes.front());
}

void remap_each(std::deque<Haplotype>& haplotypes, const GenomicRegion& region)
{
    std::transform(std::cbegin(haplotypes), std::cend(haplotypes), std::begin(haplotypes),
                   [&] (const Haplotype& haplotype) { return remap(haplotype, region); });
}

template <typename T>
void merge_unique(std::deque<T>&& src, std::deque<T>& dst)
{
    const auto itr = utils::append(std::move(src), dst);
    std::inplace_merge(std::begin(dst), itr, std::end(dst));
    dst.erase(std::unique(std::begin(dst), std::end(dst)), std::end(dst));
}

template <typename T>
void merge_unique(std::vector<T>&& src, std::vector<T>& dst)
{
    const auto itr = utils::append(std::move(src), dst);
    std::inplace_merge(std::begin(dst), itr, std::end(dst));
    dst.erase(std::unique(std::begin(dst), std::end(dst)), std::end(dst));
}

template <typename T>
auto insert_sorted(const T& val, std::deque<T>& dst)
{
    return dst.insert(std::upper_bound(std::begin(dst), std::end(dst), val), val);
}

template <typename Container>
auto find_reference(const Container& haplotypes)
{
    return std::find_if(std::cbegin(haplotypes), std::cend(haplotypes), [] (const auto& haplotype) { return is_reference(haplotype); });
}

template <typename Container>
bool has_reference(const Container& haplotypes)
{
    return find_reference(haplotypes) != std::cend(haplotypes);
}

} // namespace

std::deque<CallWrapper>
Caller::call_variants(const GenomicRegion& call_region, const MappableFlatSet<Variant>& candidates,
                      const ReadMap& reads, const ReadPipe::Report& read_report,
                      ProgressMeter& progress_meter) const
{
    std::deque<CallWrapper> result {};
    auto haplotype_likelihoods = make_haplotype_likelihood_cache();
    if (candidates.empty()) {
        if (refcalls_requested()) {
            utils::append(call_reference(call_region, reads), result);
        }
        progress_meter.log_completed(call_region);
        return result;
    }
    auto haplotype_generator = make_haplotype_generator(candidates, reads, read_report);
    GeneratorStatus status;
    std::vector<Haplotype> haplotypes {}, next_haplotypes {};
    GenomicRegion active_region;
    boost::optional<GenomicRegion> next_active_region {}, prev_called_region {}, backtrack_region {};
    auto completed_region = head_region(call_region);
    std::deque<Haplotype> protected_haplotypes {};
    while (true) {
        status = generate_active_haplotypes(call_region, haplotype_generator, active_region,
                                            next_active_region, haplotypes, next_haplotypes);
        if (status == GeneratorStatus::done) {
            if (refcalls_requested()) {
                if (!prev_called_region) {
                    utils::append(call_reference(call_region, reads), result);
                } else if (ends_before(*prev_called_region, call_region)) {
                    const auto final_refcall_region = right_overhang_region(call_region, *prev_called_region);
                    utils::append(call_reference(final_refcall_region, reads), result);
                }
            }
            progress_meter.log_completed(active_region);
            break;
        }
        if (status == GeneratorStatus::skipped) {
            haplotype_likelihoods.clear();
            progress_meter.log_completed(active_region);
            continue;
        }
        const auto active_reads = copy_overlapped(reads, active_region);
        if (!refcalls_requested() && !has_coverage(active_reads)) {
            if (debug_log_) stream(*debug_log_) << "Skipping active region " << active_region << " as there are no active reads";
            continue;
        }
        if (debug_log_) stream(*debug_log_) << "There are " << count_reads(active_reads) << " active reads in " << active_region;
        if (!populate(haplotype_likelihoods, active_region, haplotypes, candidates, active_reads)) {
            haplotype_generator.clear_progress();
            haplotype_likelihoods.clear();
            continue;
        }
        if (!protected_haplotypes.empty()) {
            assert(!haplotypes.empty());
            std::sort(std::begin(haplotypes), std::end(haplotypes));
            remap_each(protected_haplotypes, haplotype_region(haplotypes));
            std::sort(std::begin(protected_haplotypes), std::end(protected_haplotypes));
        }
        if (parameters_.protect_reference_haplotype && !has_reference(protected_haplotypes)) {
            const auto reference_haplotype_itr = find_reference(haplotypes);
            if (reference_haplotype_itr != std::cend(haplotypes)) {
                insert_sorted(*reference_haplotype_itr, protected_haplotypes);
            }
        }
        auto has_removal_impact = filter_haplotypes(haplotypes, haplotype_generator, haplotype_likelihoods, protected_haplotypes);
        if (haplotypes.empty()) continue;
        const auto caller_latents = infer_latents(haplotypes, haplotype_likelihoods);
        if (trace_log_) {
            debug::print_haplotype_posteriors(stream(*trace_log_), *caller_latents->haplotype_posteriors(), -1);
        } else if (debug_log_) {
            debug::print_haplotype_posteriors(stream(*debug_log_), *caller_latents->haplotype_posteriors());
        }
        if (!is_saturated(haplotypes, *caller_latents)) {
            filter_haplotypes(has_removal_impact, haplotypes, haplotype_generator,
                              haplotype_likelihoods, *caller_latents, protected_haplotypes);
        } else {
            if (debug_log_) {
                *debug_log_ << "Haplotypes are saturated, clearing lagging";
            }
            haplotype_generator.clear_progress();
        }
        status = generate_next_active_haplotypes(next_haplotypes, next_active_region, backtrack_region, haplotype_generator);
        if (backtrack_region) {
            // Only protect haplotypes in backtrack - or holdout - regions as these are more likely
            // to suffer from window artifacts.
            merge_unique(get_called_haplotypes(*caller_latents), protected_haplotypes);
        } else {
            protected_haplotypes.clear();
        }
        if (status != GeneratorStatus::skipped) {
            call_variants(active_region, call_region, next_active_region, backtrack_region,
                          candidates, haplotypes, haplotype_likelihoods, reads, *caller_latents,
                          result, prev_called_region, completed_region);
        }
        haplotype_likelihoods.clear();
        progress_meter.log_completed(completed_region);
    }
    return result;
}

std::size_t Caller::do_remove_duplicates(std::vector<Haplotype>& haplotypes) const
{
    return octopus::remove_duplicates(haplotypes, Haplotype {haplotype_region(haplotypes), reference_.get()});
}

boost::optional<MemoryFootprint> Caller::target_max_memory() const noexcept
{
    return parameters_.target_max_memory;
}

Caller::GeneratorStatus
Caller::generate_active_haplotypes(const GenomicRegion& call_region,
                                   HaplotypeGenerator& haplotype_generator,
                                   GenomicRegion& active_region,
                                   boost::optional<GenomicRegion>& next_active_region,
                                   std::vector<Haplotype>& haplotypes,
                                   std::vector<Haplotype>& next_haplotypes) const
{
    if (next_active_region) {
        haplotypes = std::move(next_haplotypes);
        active_region = std::move(*next_active_region);
        next_active_region = boost::none;
    } else {
        try {
            std::tie(haplotypes, next_active_region, std::ignore) = haplotype_generator.generate();
            if (next_active_region) {
                active_region = std::move(*next_active_region);
                next_active_region = boost::none;
            }
        } catch (const HaplotypeGenerator::HaplotypeOverflow& e) {
            logging::WarningLogger warn_log {};
            stream(warn_log) << "Skipping region " << e.region() << " as there are too many haplotypes";
            haplotype_generator.clear_progress();
            active_region = e.region();
            return GeneratorStatus::skipped;
        }
    }
    if (debug_log_) stream(*debug_log_) << "Active region is " << active_region;
    if (is_after(active_region, call_region) || haplotypes.empty()) {
        if (debug_log_) {
            if (haplotypes.empty()) {
                stream(*debug_log_) << "No haplotypes were generated in active region " << active_region;
            } else {
                stream(*debug_log_) << "Generated " << haplotypes.size()
                                    << " haplotypes in " << active_region
                                    << " but they are after the call region "
                                    << call_region;
            }
        }
        return GeneratorStatus::done;
    }
    remove_duplicates(haplotypes);
    if (debug_log_) stream(*debug_log_) << "There are " << haplotypes.size() << " unfiltered haplotypes";
    return GeneratorStatus::good;
}

void Caller::remove_duplicates(std::vector<Haplotype>& haplotypes) const
{
    const auto num_removed = do_remove_duplicates(haplotypes);
    if (debug_log_) {
        if (num_removed > 0) {
            stream(*debug_log_) << num_removed << " duplicate haplotypes were removed";
        } else {
            stream(*debug_log_) << "There are no duplicate haplotypes";
        }
    }
}

bool Caller::filter_haplotypes(std::vector<Haplotype>& haplotypes,
                               HaplotypeGenerator& haplotype_generator,
                               HaplotypeLikelihoodArray& haplotype_likelihoods,
                               const std::deque<Haplotype>& protected_haplotypes) const
{
    bool has_removal_impact {false};
    auto removed_haplotypes = filter(haplotypes, haplotype_likelihoods, protected_haplotypes);
    if (haplotypes.empty()) {
        // This can only happen if all haplotypes have equal likelihood
        haplotype_generator.clear_progress();
        haplotype_likelihoods.clear();
    } else {
        haplotypes.shrink_to_fit();
        haplotype_likelihoods.erase(removed_haplotypes);
        has_removal_impact = haplotype_generator.removal_has_impact();
        if (has_removal_impact) {
            haplotype_generator.remove(removed_haplotypes);
            haplotype_generator.collapse(haplotypes);
        } else {
            haplotype_generator.clear_progress();
        }
        removed_haplotypes.clear();
        removed_haplotypes.shrink_to_fit();
        if (debug_log_) stream(*debug_log_) << "There are " << haplotypes.size() << " haplotypes after filtering";
    }
    return has_removal_impact;
}

Caller::GeneratorStatus
Caller::generate_next_active_haplotypes(std::vector<Haplotype>& next_haplotypes,
                                        boost::optional<GenomicRegion>& next_active_region,
                                        boost::optional<GenomicRegion>& backtrack_region,
                                        HaplotypeGenerator& haplotype_generator) const
{
    try {
        std::tie(next_haplotypes, next_active_region, backtrack_region) = haplotype_generator.generate();
    } catch (const HaplotypeGenerator::HaplotypeOverflow& e) {
        logging::WarningLogger warn_log {};
        stream(warn_log) << "Skipping region " << e.region() << " as there are too many haplotypes";
        haplotype_generator.clear_progress();
        return GeneratorStatus::skipped;
    }
    return GeneratorStatus::good;
}

bool Caller::is_saturated(const std::vector<Haplotype>& haplotypes, const Latents& latents) const
{
    return haplotypes.size() == parameters_.max_haplotypes
           && count_probable_haplotypes(*latents.haplotype_posteriors()) > parameters_.max_haplotypes / 2;
}

unsigned Caller::count_probable_haplotypes(const Caller::Latents::HaplotypeProbabilityMap& haplotype_posteriors) const
{
    return std::count_if(std::cbegin(haplotype_posteriors), std::cend(haplotype_posteriors),
                         [this] (const auto& p) { return p.second >= parameters_.saturation_limit.probability_true(); });
}

void Caller::filter_haplotypes(bool prefilter_had_removal_impact,
                               const std::vector<Haplotype>& haplotypes,
                               HaplotypeGenerator& haplotype_generator,
                               const HaplotypeLikelihoodArray& haplotype_likelihoods,
                               const Latents& latents,
                               const std::deque<Haplotype>& protected_haplotypes) const
{
    if (prefilter_had_removal_impact) { // if there was no impact before then there can't be now either
        prefilter_had_removal_impact = haplotype_generator.removal_has_impact();
    }
    if (prefilter_had_removal_impact) {
        const auto max_to_remove = haplotype_generator.max_removal_impact();
        auto removable_haplotypes = get_removable_haplotypes(haplotypes, haplotype_likelihoods,
                                                             *latents.haplotype_posteriors(),
                                                             protected_haplotypes, max_to_remove);
        if (debug_log_) {
            stream(*debug_log_) << "Discarding " << removable_haplotypes.size()
                                << " haplotypes with low posterior support";
        }
        haplotype_generator.remove(removable_haplotypes);
    }
}

namespace {

bool have_callable_region(const GenomicRegion& active_region,
                          const boost::optional<GenomicRegion>& next_active_region,
                          const boost::optional<GenomicRegion>& backtrack_region,
                          const GenomicRegion& call_region)
{
    return (!backtrack_region || begins_before(active_region, *backtrack_region))
           && (!next_active_region || begins_before(active_region, *next_active_region))
           && overlaps(active_region, call_region);
}

auto get_passed_region(const GenomicRegion& active_region,
                       const boost::optional<GenomicRegion>& next_active_region,
                       const boost::optional<GenomicRegion>& backtrack_region)
{
    auto result = active_region;
    if (next_active_region) {
        result = *overlapped_region(left_overhang_region(result, *next_active_region), result);
    }
    if (backtrack_region) {
        result = *overlapped_region(left_overhang_region(result, *backtrack_region), result);
    }
    return result;
}

auto get_uncalled_region(const GenomicRegion& active_region,
                         const GenomicRegion& passed_region,
                         const GenomicRegion& completed_region)
{
    auto result = *overlapped_region(active_region, passed_region);
    return *overlapped_region(result, right_overhang_region(result, completed_region));
}

bool is_lhs_boundry_insertion(const Variant& variant, const GenomicRegion& region)
{
    return is_empty_region(variant) && begins_equal(variant, region);
}

bool has_noninteracting_insertion(const GenomicRegion& insertion_region,
                                  const MappableFlatSet<Variant>& candidates)
{
    auto overlapped = overlap_range(candidates, insertion_region);
    while (!empty(overlapped) && begins_before(overlapped.front(), insertion_region)) {
        overlapped.advance_begin(1);
    }
    return !empty(overlapped) && is_empty_region(overlapped.front()) && is_empty_region(overlapped.back());
}

bool can_remove_lhs_boundary_insertions(const GenomicRegion& uncalled_active_region,
                                        const boost::optional<GenomicRegion>& prev_called_region,
                                        const MappableFlatSet<Variant>& candidates)
{
    return prev_called_region && are_adjacent(*prev_called_region, uncalled_active_region)
           && has_noninteracting_insertion(head_region(uncalled_active_region), candidates);
}

bool is_rhs_boundry_insertion(const Variant& variant, const GenomicRegion& region)
{
    return is_empty_region(variant) && ends_equal(variant, region);
}

bool has_interacting_insertion(const GenomicRegion& insertion_region,
                               const MappableFlatSet<Variant>& candidates)
{
    auto overlapped = overlap_range(candidates, insertion_region);
    // Overlapped non-insertions to the left of insertion_region can be safely ignored
    // as active region boundary alleles are guaranteed not to have interacting callable alleles
    // on *both* sides. However, this does not mean there are no such candidates; there may be
    // candidate alleles on the lhs that have been removed via haplotype reduction. These must be
    // ignored.
    while (!empty(overlapped) && begins_before(overlapped.front(), insertion_region)) {
        overlapped.advance_begin(1);
    }
    return !empty(overlapped) && is_empty_region(overlapped.front()) && !is_empty_region(overlapped.back());
}

bool can_remove_rhs_boundary_insertions(const GenomicRegion& uncalled_active_region,
                                        const boost::optional<GenomicRegion>& next_active_region,
                                        const boost::optional<GenomicRegion>& backtrack_region,
                                        const MappableFlatSet<Variant>& candidates)
{
    return ((next_active_region && are_adjacent(uncalled_active_region, *next_active_region))
            || (backtrack_region && are_adjacent(uncalled_active_region, *backtrack_region)))
           && has_interacting_insertion(tail_region(uncalled_active_region), candidates);
}

template <typename R>
void remove_duplicate_boundary_insertions(R& contained, const MappableFlatSet<Variant>& candidates,
                                          const GenomicRegion& uncalled_active_region,
                                          const boost::optional<GenomicRegion>& prev_called_region,
                                          const boost::optional<GenomicRegion>& next_active_region,
                                          const boost::optional<GenomicRegion>& backtrack_region)
{
    // We only want to call insertions once, which means removing insertions that sit on the boundary
    // of an uncalled active region. Either we remove those that are at the front, or at the back.
    // It is usually better to keep those at the back as they have better posterior resolution, but
    // if there are other (non insertion) candidates overlapping with the insertion, then we must
    // keep front insertions as they must be called in phase with the interacting non-insertion.
    if (can_remove_lhs_boundary_insertions(uncalled_active_region, prev_called_region, candidates)) {
        while (!empty(contained) && is_lhs_boundry_insertion(contained.front(), uncalled_active_region)) {
            contained.advance_begin(1);
        }
    }
    if (can_remove_rhs_boundary_insertions(uncalled_active_region, next_active_region, backtrack_region, candidates)) {
        while (!empty(contained) && is_rhs_boundry_insertion(contained.back(), uncalled_active_region)) {
            contained.advance_end(-1);
        }
    }
}

auto extract_callable_variants(const MappableFlatSet<Variant>& candidates,
                               const GenomicRegion& uncalled_active_region,
                               const boost::optional<GenomicRegion>& prev_called_region,
                               const boost::optional<GenomicRegion>& next_active_region,
                               const boost::optional<GenomicRegion>& backtrack_region)
{
    auto contained = contained_range(candidates, uncalled_active_region);
    remove_duplicate_boundary_insertions(contained, candidates, uncalled_active_region, prev_called_region,
                                         next_active_region, backtrack_region);
    return std::vector<Variant> {std::cbegin(contained), std::cend(contained)};
}

} // namespace

void Caller::call_variants(const GenomicRegion& active_region,
                           const GenomicRegion& call_region,
                           const boost::optional<GenomicRegion>& next_active_region,
                           const boost::optional<GenomicRegion>& backtrack_region,
                           const MappableFlatSet<Variant>& candidates,
                           const std::vector<Haplotype>& haplotypes,
                           const HaplotypeLikelihoodArray& haplotype_likelihoods,
                           const ReadMap& reads,
                           const Latents& latents,
                           std::deque<CallWrapper>& result,
                           boost::optional<GenomicRegion>& prev_called_region,
                           GenomicRegion& completed_region) const
{
    if (have_callable_region(active_region, next_active_region, backtrack_region, call_region)) {
        const auto passed_region = get_passed_region(active_region, next_active_region, backtrack_region);
        const auto uncalled_region = get_uncalled_region(active_region, passed_region, completed_region);
        auto active_candidates = extract_callable_variants(candidates, uncalled_region, prev_called_region,
                                                           next_active_region, backtrack_region);
        std::vector<CallWrapper> calls {};
        if (!active_candidates.empty()) {
            if (debug_log_) stream(*debug_log_) << "Calling variants in region " << uncalled_region;
            calls = wrap(call_variants(active_candidates, latents));
            if (!calls.empty()) {
                set_model_posteriors(calls, latents, haplotypes, haplotype_likelihoods);
                set_phasing(calls, latents, haplotypes, call_region);
            }
        }
        if (refcalls_requested()) {
            const auto refcall_region = right_overhang_region(uncalled_region, completed_region);
            const auto pileups = make_pileups(reads, latents, refcall_region);
            auto alleles = generate_reference_alleles(refcall_region, active_candidates, calls);
            auto reference_calls = call_reference_helper(alleles, latents, pileups);
            const auto itr = utils::append(std::move(reference_calls), calls);
            std::inplace_merge(std::begin(calls), itr, std::end(calls));
        }
        utils::append(std::move(calls), result);
        prev_called_region = uncalled_region;
        completed_region = encompassing_region(completed_region, passed_region);
    }
}

Genotype<Haplotype> Caller::call_genotype(const Latents& latents, const SampleName& sample) const
{
    const auto genotype_posteriors_ptr = latents.genotype_posteriors();
    const auto& genotype_posteriors = (*genotype_posteriors_ptr)[sample];
    auto itr = std::max_element(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                                [] (const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
    assert(itr != std::cend(genotype_posteriors));
    return itr->first;
}

std::deque<Haplotype> Caller::get_called_haplotypes(const Latents& latents) const
{
    std::deque<Haplotype> result {};
    for (const auto& sample : samples_) {
        auto called_genotype = call_genotype(latents, sample);
        std::copy(std::begin(called_genotype), std::cend(called_genotype), std::back_inserter(result));
    }
    std::sort(std::begin(result), std::end(result));
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

bool requires_model_evaluation(const std::vector<CallWrapper>& calls)
{
    return std::any_of(std::cbegin(calls), std::cend(calls),
                       [] (const auto& call) { return call->requires_model_evaluation(); });
}

void Caller::set_model_posteriors(std::vector<CallWrapper>& calls, const Latents& latents,
                                  const std::vector<Haplotype>& haplotypes,
                                  const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    if (parameters_.allow_model_filtering || requires_model_evaluation(calls)) {
        const auto mp = calculate_model_posterior(haplotypes, haplotype_likelihoods, latents);
        if (mp) {
            for (auto& call : calls) {
                call->set_model_posterior(probability_to_phred(1 - *mp));
            }
        }
    }
}

namespace {

auto get_phase_regions(const MappableFlatSet<Variant>& candidates,
                       const GenomicRegion& active_region)
{
    return extract_regions(contained_range(candidates, active_region));
}

void set_phase(const SampleName& sample, const Phaser::PhaseSet::PhaseRegion& phase,
               const std::vector<GenomicRegion>& call_regions, CallWrapper& call)
{
    const auto overlapped = overlap_range(call_regions, phase.mapped_region(), BidirectionallySortedTag {});
    if (!overlapped.empty()) {
        call->set_phase(sample, Call::PhaseCall {encompassing_region(overlapped.front(), call), phase.score});
    } else {
        call->set_phase(sample, Call::PhaseCall {mapped_region(call), phase.score});
    }
}

void set_phasing(std::vector<CallWrapper>& calls, const Phaser::PhaseSet& phase_set, const GenomicRegion& calling_region)
{
    if (!calls.empty()) {
        const auto call_regions = extract_regions(calls);
        for (auto& call : calls) {
            const auto& call_region = mapped_region(call);
            for (const auto& p : phase_set.phase_regions) {
                const SampleName& sample {p.first};
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
                            set_phase(sample, clipped_phase, call_regions, call);
                        }
                    } else {
                        set_phase(sample, *phase, call_regions, call);
                    }
                }
            }
        }
    }
}

} // namespace

void Caller::set_phasing(std::vector<CallWrapper>& calls, const Latents& latents,
                         const std::vector<Haplotype>& haplotypes,
                         const GenomicRegion& call_region) const
{
    const auto phase = phaser_.force_phase(haplotypes, *latents.genotype_posteriors(),
                                           extract_regions(calls), get_genotype_calls(latents));
    if (debug_log_) debug::print_phase_sets(stream(*debug_log_), phase);
    octopus::set_phasing(calls, phase, call_region);
}

bool Caller::refcalls_requested() const noexcept
{
    return parameters_.refcall_type != RefCallType::none;
}

bool check_reference(const Variant& v, const ReferenceGenome& reference)
{
    return ref_sequence(v) == reference.fetch_sequence(mapped_region(v));
}

bool check_reference(const std::vector<Variant>& variants, const ReferenceGenome& reference)
{
    return std::all_of(std::cbegin(variants), std::cend(variants), [&] (const auto& v) { return check_reference(v, reference); });
}

MappableFlatSet<Variant> Caller::generate_candidate_variants(const GenomicRegion& region) const
{
    if (debug_log_) stream(*debug_log_) << "Generating candidate variants in region " << region;
    auto raw_candidates = candidate_generator_.generate(region);
    if (debug_log_) debug::print_left_aligned_candidates(stream(*debug_log_), raw_candidates, reference_);
    auto final_candidates = unique_left_align(std::move(raw_candidates), reference_);
    assert(check_reference(final_candidates, reference_));
    candidate_generator_.clear();
    return MappableFlatSet<Variant> {std::make_move_iterator(std::begin(final_candidates)),
                                     std::make_move_iterator(std::end(final_candidates))};
}

HaplotypeGenerator Caller::make_haplotype_generator(const MappableFlatSet<Variant>& candidates,
                                                    const ReadMap& reads, const ReadPipe::Report& read_report) const
{
    return haplotype_generator_builder_.build(reference_, candidates, reads, read_report);
}

HaplotypeLikelihoodArray Caller::make_haplotype_likelihood_cache() const
{
    return HaplotypeLikelihoodArray {likelihood_model_, parameters_.max_haplotypes, samples_};
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
        lhs_flank = expand_rhs(lhs_flank, -1); // stops boundary insertions being inactive
    }
    const auto lhs_inactive_candidates = contained_range(candidates, lhs_flank);
    if (lhs_inactive_candidates.empty()) {
        lhs_flank = head_region(lhs_flank);
    } else {
        lhs_flank = closed_region(lhs_flank, rightmost_region(lhs_inactive_candidates));
    }
    if (is_empty_region(*rightmost_mappable(active_candidates)) && !is_empty(rhs_flank)) {
        rhs_flank = expand_lhs(rhs_flank, -1); // stops boundary insertions being inactive
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
    const auto flank_regions = calculate_flank_regions(haplotype_region(haplotypes), active_region, candidates);
    return HaplotypeLikelihoodArray::FlankState {
        size(flank_regions.first), size(flank_regions.second)
    };
}

bool Caller::populate(HaplotypeLikelihoodArray& haplotype_likelihoods,
                      const GenomicRegion& active_region,
                      const std::vector<Haplotype>& haplotypes,
                      const MappableFlatSet<Variant>& candidates,
                      const ReadMap& active_reads) const
{
    assert(haplotype_likelihoods.is_empty());
    boost::optional<HaplotypeLikelihoodArray::FlankState> flank_state {};
    if (debug_log_) {
        stream(*debug_log_) << "Calculating likelihoods for " << haplotypes.size() << " haplotypes";
        debug::print_active_candidates(stream(*debug_log_), candidates, active_region);
        stream(*debug_log_) << "Haplotype region is " << haplotype_region(haplotypes);
    }
    if (likelihood_model_.can_use_flank_state()) {
        flank_state = calculate_flank_state(haplotypes, active_region, candidates);
        if (debug_log_) {
            debug::print_inactive_flanking_candidates(stream(*debug_log_), candidates, active_region,
                                                      haplotype_region(haplotypes));
        }
    }
    try {
        haplotype_likelihoods.populate(active_reads, haplotypes, std::move(flank_state));
    } catch(const HaplotypeLikelihoodModel::ShortHaplotypeError& e) {
        if (debug_log_) {
            stream(*debug_log_) << "Skipping " << active_region << " as a haplotype was too short by "
                                << e.required_extension() << "bp";
        }
        return false;
    }
    if (trace_log_) {
        debug::print_read_haplotype_likelihoods(stream(*trace_log_), haplotypes, active_reads,
                                               haplotype_likelihoods, -1);
    }
    return true;
}

std::vector<Haplotype>
Caller::filter(std::vector<Haplotype>& haplotypes, const HaplotypeLikelihoodArray& haplotype_likelihoods,
               const std::deque<Haplotype>& protected_haplotypes) const
{
    std::vector<Haplotype> removed_haplotypes {};
    if (protected_haplotypes.empty()) {
        removed_haplotypes = filter_to_n(haplotypes, samples_, haplotype_likelihoods, parameters_.max_haplotypes);
    } else {
        if (debug_log_) {
            stream(*debug_log_) << "Protecting " << protected_haplotypes.size() << " haplotypes from filtering";
        }
        std::vector<Haplotype> removable_haplotypes {}, protected_copies {};
        removable_haplotypes.reserve(haplotypes.size());
        protected_copies.reserve(protected_haplotypes.size());
        assert(std::is_sorted(std::cbegin(haplotypes), std::cend(haplotypes)));
        assert(std::is_sorted(std::cbegin(protected_haplotypes), std::cend(protected_haplotypes)));
        std::set_difference(std::cbegin(haplotypes), std::cend(haplotypes),
                            std::cbegin(protected_haplotypes), std::cend(protected_haplotypes),
                            std::back_inserter(removable_haplotypes));
        std::set_intersection(std::cbegin(haplotypes), std::cend(haplotypes),
                              std::cbegin(protected_haplotypes), std::cend(protected_haplotypes),
                              std::back_inserter(protected_copies));
        removed_haplotypes = filter_to_n(removable_haplotypes, samples_, haplotype_likelihoods, parameters_.max_haplotypes);
        haplotypes = std::move(removable_haplotypes);
        std::sort(std::begin(haplotypes), std::end(haplotypes));
        merge_unique(std::move(protected_copies), haplotypes);
    }
    const auto reference_itr = std::find_if(std::begin(removed_haplotypes), std::end(removed_haplotypes),
                                            [] (const auto& haplotype) { return is_reference(haplotype); });
    if (reference_itr != std::end(removed_haplotypes)) {
        if (debug_log_) {
            stream(*debug_log_) << "Putting back filtered reference haplotype";
        }
        haplotypes.push_back(std::move(*reference_itr));
        removed_haplotypes.erase(reference_itr);
    }
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
                                 const HaplotypeLikelihoodArray& haplotype_likelihoods,
                                 const Caller::Latents::HaplotypeProbabilityMap& haplotype_posteriors,
                                 const std::deque<Haplotype>& protected_haplotypes, const unsigned max_to_remove) const
{
    if (debug_log_) {
        stream(*debug_log_) << "Protecting " << protected_haplotypes.size() << " haplotypes from filtering";
    }
    if (protected_haplotypes.empty()) {
        return extract_removable(haplotypes, haplotype_posteriors, samples_, haplotype_likelihoods,
                                 max_to_remove, parameters_.haplotype_extension_threshold.probability_false());
    } else if (protected_haplotypes.size() == 1 && parameters_.protect_reference_haplotype && is_reference(protected_haplotypes.front())) {
        auto result = extract_removable(haplotypes, haplotype_posteriors, samples_, haplotype_likelihoods,
                                        max_to_remove, parameters_.haplotype_extension_threshold.probability_false());
        auto reference_itr = find_reference(result);
        if (reference_itr != std::cend(result)) result.erase(reference_itr);
        return result;
    } else {
        std::vector<Haplotype> removable_haplotypes {};
        removable_haplotypes.reserve(haplotypes.size());
        assert(std::is_sorted(std::cbegin(haplotypes), std::cend(haplotypes)));
        assert(std::is_sorted(std::cbegin(protected_haplotypes), std::cend(protected_haplotypes)));
        std::set_difference(std::cbegin(haplotypes), std::cend(haplotypes),
                            std::cbegin(protected_haplotypes), std::cend(protected_haplotypes),
                            std::back_inserter(removable_haplotypes));
        if (debug_log_) {
            stream(*debug_log_) << "There are " << removable_haplotypes.size() << " removable haplotypes";
        }
        return extract_removable(removable_haplotypes, haplotype_posteriors, samples_, haplotype_likelihoods,
                                 max_to_remove, parameters_.haplotype_extension_threshold.probability_false());
    }
}

Caller::GenotypeCallMap Caller::get_genotype_calls(const Latents& latents) const
{
    GenotypeCallMap result {samples_.size()};
    for (const auto& sample : samples_) {
        result.emplace(sample, call_genotype(latents, sample));
    }
    return result;
}

bool Caller::done_calling(const GenomicRegion& region) const noexcept
{
    return is_empty(region);
}

bool Caller::is_merge_block_refcalling() const noexcept
{
    return parameters_.refcall_type == RefCallType::blocked && parameters_.refcall_block_merge_threshold;
}

std::vector<CallWrapper> Caller::call_reference(const GenomicRegion& region, const ReadMap& reads) const
{
    const auto active_reads = copy_overlapped(reads, region);
    auto haplotype_likelihoods = make_haplotype_likelihood_cache();
    std::vector<Haplotype> haplotypes;
    if (has_coverage(active_reads)) {
        const auto active_reads_region = encompassing_region(active_reads);
        const auto haplotype_region = expand(active_reads_region, HaplotypeLikelihoodModel{}.pad_requirement());
        haplotypes.emplace_back(haplotype_region, reference_);
    } else {
        haplotypes.emplace_back(region, reference_);
    }
    haplotype_likelihoods.populate(active_reads, haplotypes);
    const auto latents = infer_latents(haplotypes, haplotype_likelihoods);
    const auto pileups = make_pileups(active_reads, *latents, region);
    const auto alleles = generate_reference_alleles(region);
    return call_reference_helper(alleles, *latents, pileups);
}

std::vector<CallWrapper>
Caller::call_reference_helper(const std::vector<Allele>& alleles, const Latents& latents, const ReadPileupMap& pileups) const
{
    if (is_merge_block_refcalling()) {
        return wrap(squash_reference_calls(call_reference(alleles, latents, pileups)));
    } else {
        return wrap(call_reference(alleles, latents, pileups));
    }
}

namespace {

template <typename Container>
auto extract_unique_regions(const Container& mappables)
{
    auto result = extract_regions(mappables);
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

template <typename T1, typename T2>
auto set_difference(std::vector<T1>&& first, const std::vector<T2>& second)
{
    std::vector<T1> result {};
    result.reserve(first.size());
    std::set_difference(std::make_move_iterator(std::begin(first)),
                        std::make_move_iterator(std::end(first)),
                        std::cbegin(second), std::cend(second),
                        std::back_inserter(result));
    return result;
}

auto extract_uncalled_candidate_regions(const std::vector<Variant>& candidates,
                                        const std::vector<CallWrapper>& calls)
{
    auto uncalled_regions = set_difference(extract_unique_regions(candidates), extract_unique_regions(calls));
    return extract_covered_regions(std::move(uncalled_regions));
}

template <typename T1, typename T2>
auto merge(std::vector<T1>&& first, std::vector<T2>&& second)
{
    std::vector<GenomicRegion> result {};
    result.reserve(first.size() + second.size());
    using std::make_move_iterator; using std::begin; using std::end;
    std::merge(make_move_iterator(begin(first)), make_move_iterator(end(first)),
               make_move_iterator(begin(second)), make_move_iterator(end(second)),
               std::back_inserter(result));
    return result;
}

auto extract_uncalled_reference_regions(const GenomicRegion& region,
                                        const std::vector<Variant>& candidates,
                                        const std::vector<CallWrapper>& calls)
{
    auto uncalled_candidate_regions = extract_uncalled_candidate_regions(candidates, calls);
    auto noncandidate_regions = extract_intervening_regions(extract_covered_regions(candidates), region);
    auto result = merge(std::move(uncalled_candidate_regions), std::move(noncandidate_regions));
    result.erase(std::remove_if(std::begin(result), std::end(result),
                                [] (const auto& region) { return is_empty(region); }),
                 std::end(result));
    return result;
}

auto make_positional_reference_alleles(const std::vector<GenomicRegion>& regions, const ReferenceGenome& reference)
{
    std::vector<Allele> result {};
    result.reserve(sum_region_sizes(regions));
    for (const auto& region : regions) {
        utils::append(make_positional_reference_alleles(region, reference), result);
    }
    return result;
}

} // namespace

std::vector<Allele>
Caller::generate_reference_alleles(const GenomicRegion& region,
                                   const std::vector<Variant>& candidates,
                                   const std::vector<CallWrapper>& calls) const
{
    auto refcall_regions = extract_uncalled_reference_regions(region, candidates, calls);
    if (parameters_.refcall_type == RefCallType::positional || is_merge_block_refcalling()) {
        return make_positional_reference_alleles(std::move(refcall_regions), reference_);
    } else {
        return make_reference_alleles(std::move(refcall_regions), reference_);
    }
}

std::vector<Allele> Caller::generate_reference_alleles(const GenomicRegion& region) const
{
    return generate_reference_alleles(region, {}, {});
}

namespace {

auto overlap_range(std::vector<ReadPileup>& pileups, const AlignedRead& read)
{
    return overlap_range(std::begin(pileups), std::end(pileups), contig_region(read), BidirectionallySortedTag {});
}

} // namespace

auto make_pileups(const std::vector<AlignedRead>& reads, const Genotype<Haplotype>& genotype, const GenomicRegion& region)
{
    const auto realignments = assign_and_realign(reads, genotype);
    ReadPileups result {};
    result.reserve(size(region));
    for (auto position = region.begin(); position < region.end(); ++position) {
        result.emplace_back(position);
    }
    for (const auto& p : realignments) {
        for (const auto& read : p.second) {
            for (ReadPileup& pileup : overlap_range(result, read)) {
                pileup.add(read);
            }
        }
    }
    return result;
}

ReadPileups make_pileups(const ReadContainer& reads, const Genotype<Haplotype>& genotype, const GenomicRegion& region)
{
    const auto overlapped_reads = overlap_range(reads, region);
    const std::vector<AlignedRead> active_reads {std::cbegin(overlapped_reads), std::cend(overlapped_reads)};
    if (!active_reads.empty()) {
        const auto active_reads_region = encompassing_region(active_reads);
        const auto min_genotype_region = expand(active_reads_region, max_read_length(active_reads));
        if (contains(genotype, min_genotype_region)) {
            return make_pileups(active_reads, genotype, region);
        } else {
            const auto expanded_genotype = remap(genotype, min_genotype_region);
            return make_pileups(active_reads, expanded_genotype, region);
        }
    } else {
        return make_pileups(active_reads, genotype, region);
    }
}

Caller::ReadPileupMap Caller::make_pileups(const ReadMap& reads, const Latents& latents, const GenomicRegion& region) const
{
    ReadPileupMap result {};
    result.reserve(samples_.size());
    for (const auto& sample : samples_) {
        const auto called_genotype = call_genotype(latents, sample);
        result.emplace(sample, octopus::make_pileups(reads.at(sample), called_genotype, region));
    }
    return result;
}

namespace {

auto get_min_quality(const std::vector<std::unique_ptr<ReferenceCall>>& refcalls)
{
    auto itr = std::min_element(std::cbegin(refcalls), std::cend(refcalls),
                                [] (const auto& lhs, const auto& rhs) { return lhs->quality() < rhs->quality(); });
    return (*itr)->quality();
}

std::unique_ptr<ReferenceCall>
concat(const std::vector<std::unique_ptr<ReferenceCall>>& refcalls, const std::vector<SampleName>& samples)
{
    assert(!refcalls.empty());
    auto region = encompassing_region(refcalls.front()->mapped_region(), refcalls.back()->mapped_region());
    Allele::NucleotideSequence sequence {};
    sequence.reserve(size(region));
    for (const auto& refcall : refcalls) {
        sequence.insert(std::cend(sequence), std::cbegin(refcall->reference().sequence()), std::cend(refcall->reference().sequence()));
    }
    Allele reference {std::move(region), std::move(sequence)};
    const auto quality = get_min_quality(refcalls);
    std::map<SampleName, ReferenceCall::GenotypeCall> genotypes {};
    for (const auto& sample : samples) {
        const auto& genotype_call = refcalls.front()->get_genotype_call(sample);
        genotypes.emplace(sample, ReferenceCall::GenotypeCall {genotype_call.genotype.ploidy(), genotype_call.posterior});
    }
    return std::make_unique<ReferenceCall>(std::move(reference), quality, std::move(genotypes));
}

bool are_similar_quality(const ReferenceCall& lhs, const ReferenceCall& rhs, const Phred<double> threshold) noexcept
{
    return std::abs(lhs.quality().score() - rhs.quality().score()) < threshold.score();
}

} // namespace

std::vector<std::unique_ptr<ReferenceCall>>
Caller::squash_reference_calls(std::vector<std::unique_ptr<ReferenceCall>> refcalls) const
{
    assert(parameters_.refcall_block_merge_threshold);
    std::vector<std::unique_ptr<ReferenceCall>> result {}, buffer {};
    if (refcalls.empty()) return result;
    result.reserve(refcalls.size());
    buffer.reserve(refcalls.size());
    for (auto& refcall : refcalls) {
        if (buffer.empty() || (are_adjacent(*buffer.back(), *refcall)
                && are_similar_quality(*buffer.front(), *refcall, *parameters_.refcall_block_merge_threshold))) {
            buffer.push_back(std::move(refcall));
        } else {
            result.push_back(concat(buffer, samples_));
            buffer.clear();
            buffer.push_back(std::move(refcall));
        }
    }
    result.push_back(concat(buffer, samples_));
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
void print_final_candidates(S&& stream, const MappableFlatSet<Variant>& candidates, const GenomicRegion& region,
                            bool number_only)
{
    if (candidates.empty()) {
        stream << "There are no final candidates in " << region << '\n';
    } else {
        if (candidates.size() == 1) {
            stream << "There is 1 final candidates in " << region << ":" << '\n';
        } else {
            stream << "There are " << candidates.size() << " final candidates in " << region << ":" << '\n';
        }
        if (!number_only) {
            for (const auto& c : candidates) stream << c << '\n';
        }
    }
}

void print_final_candidates(const MappableFlatSet<Variant>& candidates, const GenomicRegion& region, bool number_only)
{
    print_final_candidates(std::cout, candidates, region, number_only);
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
    HaplotypeLikelihoodArray cache {1, {test_sample}};
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
//    auto haplotype = make_haplotype(haplotype_str, haplotype_region_str, reference);
//    std::cout << "Haplotype: " << haplotype << std::endl;
//    debug::print_variant_alleles(haplotype);
//    std::cout << std::endl;
//    const auto active_region = io::parse_region(active_region_str, reference);
//    auto flank_state = calculate_flank_state({haplotype}, active_region, candidates);
//    std::cout << "Flank sizes: " << flank_state.lhs_flank << " " << flank_state.rhs_flank << std::endl;
//    auto read = *find_first_read(read_region_str, cigar_str, reads);
//    std::cout << "Read: " << read << std::endl;
//    auto likelihood = calculate_likelihood(haplotype, read, flank_state);
//    std::cout << "Likelihood = " << likelihood << std::endl;
}

} // namespace debug

} // namespace octopus
