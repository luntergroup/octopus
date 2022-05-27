// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "caller.hpp"

#include <algorithm>
#include <utility>
#include <tuple>
#include <iterator>
#include <stdexcept>
#include <cassert>
#include <iostream>
#include <limits>
#include <queue>

#include "concepts/mappable.hpp"
#include "basics/aligned_template.hpp"
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
#include "utils/erase_if.hpp"
#include "utils/map_utils.hpp"
#include "utils/thread_pool.hpp"

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
, bad_region_detector_ {std::move(components.bad_region_detector)}
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

std::deque<VcfRecord> 
Caller::call(const GenomicRegion& call_region, 
             ProgressMeter& progress_meter,
             OptionalThreadPool workers) const
{
    ReadPipe::Report reads_report {};
    ReadMap reads;
    boost::optional<TemplateMap> read_templates {};
    if (candidate_generator_.requires_reads()) {
        reads = read_pipe_.get().fetch_reads(expand(call_region, 100), reads_report);
        read_templates = make_read_templates(reads);
        if (read_templates) {
            add_reads(*read_templates, candidate_generator_);
        } else {
            add_reads(reads, candidate_generator_);
        }
        if (!refcalls_requested() && all_empty(reads)) {
            if (debug_log_) stream(*debug_log_) << "Stopping early as no reads found in call region " << call_region;
            return {};
        }
        if (debug_log_) stream(*debug_log_) << "Using " << count_reads(reads) << " reads in call region " << call_region;
    }
    const auto candidate_region = calculate_candidate_region(call_region, reads, reference_, candidate_generator_);
    auto candidates = generate_candidate_variants(candidate_region, workers);
    if (debug_log_) debug::print_final_candidates(stream(*debug_log_), candidates, candidate_region);
    if (!refcalls_requested() && candidates.empty()) {
        progress_meter.log_completed(call_region);
        return {};
    }
    if (!candidate_generator_.requires_reads()) {
        // as we didn't fetch them earlier
        reads = read_pipe_.get().fetch_reads(call_region, reads_report);
        read_templates = make_read_templates(reads);
    }
    std::vector<GenomicRegion> likely_difficult_regions {};
    if (bad_region_detector_ && has_coverage(reads)) {
        const auto bad_regions = bad_region_detector_->detect(candidates, reads, reads_report);
        for (const auto& bad_region : bad_regions) {
            if (bad_region.severity == BadRegionDetector::BadRegion::Severity::high) {
                if (debug_log_) {
                    stream(*debug_log_) << "Erasing " << count_contained(candidates, bad_region.region)
                                        << " candidate variants in bad region " << bad_region.region;
                }
                candidates.erase_contained(bad_region.region);
            } else{
                likely_difficult_regions.push_back(bad_region.region);
            }
        }
        likely_difficult_regions.shrink_to_fit();
    }
    auto haplotype_generator = make_haplotype_generator(candidates, reads, read_templates);
    for (auto& region : likely_difficult_regions) haplotype_generator.add_lagging_exclusion_zone(region);
    auto calls = call_variants(call_region, candidates, reads, read_templates, haplotype_generator, progress_meter, workers);
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
    for (auto& p : result) {
        realign_to_reference(p.second, p.first);
        std::sort(std::begin(p.second), std::end(p.second));
    }
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

template <typename S, typename Container>
void print_haplotypes(S&& stream, const Container& haplotypes,
                      Resolution resolution = Resolution::sequenceAndAlleles);
template <typename Container>
void print_haplotypes(const Container& haplotypes,
                      Resolution resolution = Resolution::sequenceAndAlleles);

template <typename S, typename Map>
void print_haplotype_posteriors(S&& stream, const Map& haplotype_posteriors, std::size_t n = std::numeric_limits<std::size_t>::max());
template <typename Map>
void print_haplotype_posteriors(const Map& haplotype_posteriors, std::size_t n = std::numeric_limits<std::size_t>::max());

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

template <typename Container1, typename MappableType, typename Container2>
auto append(Container1&& src, MappableBlock<MappableType, Container2>& dst)
{
    return utils::append(std::forward<Container1>(src), static_cast<Container2&>(dst));
}

namespace {

void remap_each(std::deque<Haplotype>& haplotypes, const GenomicRegion& region)
{
    std::transform(std::cbegin(haplotypes), std::cend(haplotypes), std::begin(haplotypes),
                   [&] (const Haplotype& haplotype) { return remap(haplotype, region); });
}

template <typename Container1, typename Container2>
void merge_unique(Container1&& src, Container2& dst)
{
    using utils::append;
    const auto itr = append(std::move(src), dst);
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

TemplateContainer make_read_templates_helper(const ReadContainer& reads, const ReadLinkageConfig& linkage)
{
    std::vector<AlignedTemplate> buffer {};
    buffer.reserve(reads.size());
    make_read_templates(std::cbegin(reads), std::cend(reads), std::back_inserter(buffer), linkage);
    std::sort(std::begin(buffer), std::end(buffer));
    return {ForwardSortedTag {}, std::make_move_iterator(std::begin(buffer)), std::make_move_iterator(std::end(buffer))};
}

TemplateMap make_read_templates_helper(const ReadMap& reads, const ReadLinkageConfig& linkage)
{
    TemplateMap result {};
    result.reserve(reads.size());
    for (const auto& p : reads) {
        result.emplace(p.first, make_read_templates_helper(p.second, linkage));
    }
    return result;
}

AlignedTemplate copy_overlapped(const AlignedTemplate& reads, const GenomicRegion& region)
{
    const auto overlapped = overlap_range(reads, region);
    return AlignedTemplate {std::cbegin(overlapped), std::cend(overlapped)};
}

TemplateContainer copy_overlapped(const TemplateContainer& templates, const GenomicRegion& region)
{
    const auto overlapped = overlap_range(templates, region);
    TemplateContainer result {};
    result.reserve(size(overlapped));
    for (const auto& read_template : overlapped) {
        // Need to check overlap of reads in template as template region is encompassing
        // region of all read fragements
        if (has_overlapped(read_template, region)) {
            result.insert(copy_overlapped(read_template, region));
        }
    }
    return result;
}

TemplateMap copy_overlapped(const TemplateMap& templates, const GenomicRegion& region)
{
    TemplateMap result {};
    result.reserve(templates.size());
    for (const auto& p : templates) {
        result.emplace(p.first, copy_overlapped(p.second, region));
    }
    return result;
}

bool has_coverage(const boost::variant<ReadMap, TemplateMap>& reads)
{
    return boost::apply_visitor([] (const auto& reads) { return octopus::has_coverage(reads); }, reads);
}

std::size_t count_reads(const boost::variant<ReadMap, TemplateMap>& reads)
{
    return boost::apply_visitor([] (const auto& reads) { return octopus::count_reads(reads); }, reads);
}

bool have_callable_region(const GenomicRegion& active_region,
                          const boost::optional<GenomicRegion>& next_active_region,
                          const boost::optional<GenomicRegion>& backtrack_region,
                          const GenomicRegion& call_region)
{
    return (!backtrack_region || begins_before(active_region, *backtrack_region))
           && (!next_active_region || begins_before(active_region, *next_active_region))
           && overlaps(active_region, call_region);
}

} // namespace

boost::optional<TemplateMap> Caller::make_read_templates(const ReadMap& reads) const
{
    if (parameters_.read_linkage != ReadLinkageType::none) {
        ReadLinkageConfig config {};
        config.linkage = parameters_.read_linkage;
        config.max_insert_size = 10'000;
        return make_read_templates_helper(reads, config);
    } else {
        return boost::none;
    }
}

std::deque<CallWrapper>
Caller::call_variants(const GenomicRegion& call_region,
                      const MappableFlatSet<Variant>& candidates,
                      const ReadMap& reads,
                      const boost::optional<TemplateMap>& read_templates,
                      HaplotypeGenerator& haplotype_generator,
                      ProgressMeter& progress_meter,
                      OptionalThreadPool workers) const
{
    auto haplotype_likelihoods = make_haplotype_likelihood_cache();
    std::deque<CallWrapper> result {};
    if (candidates.empty()) {
        if (refcalls_requested()) {
            utils::append(call_reference(call_region, reads), result);
        }
        progress_meter.log_completed(call_region);
        return result;
    }
    if (haplotype_generator.done()) {
        logging::WarningLogger warn_log {};
        stream(warn_log) << "No variants were considered in " << call_region << " as the region was considered uncallable";
        if (refcalls_requested()) {
            utils::append(call_reference(call_region, reads), result);
        }
        progress_meter.log_completed(call_region);
        return result;
    }
    GeneratorStatus status;
    HaplotypeBlock haplotypes {call_region}, next_haplotypes {call_region};
    GenomicRegion active_region;
    boost::optional<GenomicRegion> next_active_region {}, prev_called_region {}, backtrack_region {};
    auto completed_region = head_region(call_region);
    std::deque<Haplotype> protected_haplotypes {};
    boost::variant<ReadMap, TemplateMap> active_reads;
    while (true) {
        status = generate_active_haplotypes(call_region, haplotype_generator, active_region, next_active_region,
                                            haplotypes, next_haplotypes, backtrack_region);
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
        if (read_templates) {
            active_reads = copy_overlapped(*read_templates, active_region);
        } else {
            active_reads = copy_overlapped(reads, active_region);
        }
        if (!refcalls_requested() && !has_coverage(active_reads)) {
            if (debug_log_) stream(*debug_log_) << "Skipping active region " << active_region << " as there are no active reads";
            continue;
        }
        if (debug_log_) stream(*debug_log_) << "There are " << count_reads(active_reads) << " active reads in " << active_region;
        if (!compute_haplotype_likelihoods(haplotype_likelihoods, active_region, haplotypes, candidates, active_reads, workers)) {
            haplotype_generator.clear_progress();
            haplotype_likelihoods.clear();
            continue;
        }
        if (!protected_haplotypes.empty()) {
            assert(!haplotypes.empty());
            std::sort(std::begin(haplotypes), std::end(haplotypes));
            remap_each(protected_haplotypes, mapped_region(haplotypes));
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
        const auto caller_latents = infer_latents(haplotypes, haplotype_likelihoods, workers);
        if (trace_log_) {
            debug::print_haplotype_posteriors(stream(*trace_log_), *caller_latents->haplotype_posteriors());
        } else if (debug_log_) {
            debug::print_haplotype_posteriors(stream(*debug_log_), *caller_latents->haplotype_posteriors(), 10);
        }
        if (!is_saturated(haplotypes, *caller_latents)) {
            filter_haplotypes(has_removal_impact, haplotypes, haplotype_generator,
                              haplotype_likelihoods, *caller_latents, protected_haplotypes);
        } else {
            if (debug_log_) *debug_log_ << "Haplotypes are saturated, clearing lagging";
            haplotype_generator.clear_progress();
        }
        if (try_early_detect_phase_regions(haplotypes, candidates, active_region, *caller_latents, backtrack_region)) {
            auto phased_region = find_phased_head(haplotypes, candidates, active_region, *caller_latents);
            if (phased_region) {
                if (debug_log_) stream(*debug_log_) << "Detected phased region " << *phased_region << "... removing this from haplotype generator";
                haplotype_generator.remove(*phased_region);
            }
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
            if (have_callable_region(active_region, next_active_region, backtrack_region, call_region)) {
                call_variants(active_region, call_region, next_active_region, backtrack_region,
                              candidates, haplotypes, haplotype_likelihoods, reads, *caller_latents,
                              result, prev_called_region, completed_region);
            }
        }
        haplotype_likelihoods.clear();
        progress_meter.log_completed(completed_region);
    }
    return result;
}

std::size_t Caller::do_remove_duplicates(HaplotypeBlock& haplotypes) const
{
    return octopus::remove_duplicates(haplotypes, Haplotype {mapped_region(haplotypes), reference_.get()});
}

boost::optional<MemoryFootprint> Caller::target_max_memory() const noexcept
{
    return parameters_.target_max_memory;
}

ExecutionPolicy Caller::exucution_policy() const noexcept
{
    return parameters_.execution_policy;
}

Caller::GeneratorStatus
Caller::generate_active_haplotypes(const GenomicRegion& call_region,
                                   HaplotypeGenerator& haplotype_generator,
                                   GenomicRegion& active_region,
                                   boost::optional<GenomicRegion>& next_active_region,
                                   HaplotypeBlock& haplotypes,
                                   HaplotypeBlock& next_haplotypes,
                                   boost::optional<GenomicRegion> backtrack_region) const
{
    if (next_active_region) {
        haplotypes = std::move(next_haplotypes);
        active_region = std::move(*next_active_region);
        next_active_region = boost::none;
    } else {
        try {
            auto packet = haplotype_generator.generate();
            haplotypes = std::move(packet.haplotypes);
            if (packet.active_region) {
                active_region = std::move(*packet.active_region);
            }
            next_active_region = boost::none;
            backtrack_region = std::move(packet.backtrack_region);
        } catch (const HaplotypeGenerator::HaplotypeOverflow& e) {
            logging::WarningLogger warn_log {};
            stream(warn_log) << "Skipping region " << e.region() << " as there are too many haplotypes";
            haplotype_generator.clear_progress();
            active_region = e.region();
            return GeneratorStatus::skipped;
        }
    }
    if (debug_log_) stream(*debug_log_) << "Active region is " << active_region;
    if ((is_after(active_region, call_region) && !backtrack_region) || haplotypes.empty()) {
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
    if (debug_log_) stream(*debug_log_) << "Generated " << haplotypes.size() << " haplotypes in " << mapped_region(haplotypes);
    std::sort(std::begin(haplotypes), std::end(haplotypes));
    remove_duplicates(haplotypes);
    if (debug_log_) stream(*debug_log_) << "There are " << haplotypes.size() << " haplotypes after removing duplicates";
    return GeneratorStatus::good;
}

void Caller::remove_duplicates(HaplotypeBlock& haplotypes) const
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

bool Caller::filter_haplotypes(HaplotypeBlock& haplotypes,
                               HaplotypeGenerator& haplotype_generator,
                               HaplotypeLikelihoodArray& haplotype_likelihoods,
                               const std::deque<Haplotype>& protected_haplotypes) const
{
    bool has_removal_impact {false};
    auto removed_haplotypes = filter(haplotypes, haplotype_likelihoods, protected_haplotypes);
    std::sort(std::begin(haplotypes), std::end(haplotypes));
    if (haplotypes.empty()) {
        // This can only happen if all haplotypes have equal likelihood
        haplotype_generator.clear_progress();
        haplotype_likelihoods.clear();
    } else {
        haplotypes.shrink_to_fit();
        haplotype_likelihoods.reset(haplotypes);
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
Caller::generate_next_active_haplotypes(HaplotypeBlock& next_haplotypes,
                                        boost::optional<GenomicRegion>& next_active_region,
                                        boost::optional<GenomicRegion>& backtrack_region,
                                        HaplotypeGenerator& haplotype_generator) const
{
    try {
        auto packet = haplotype_generator.generate();
        next_haplotypes = std::move(packet.haplotypes);
        next_active_region = std::move(packet.active_region);
        backtrack_region = std::move(packet.backtrack_region);
    } catch (const HaplotypeGenerator::HaplotypeOverflow& e) {
        logging::WarningLogger warn_log {};
        stream(warn_log) << "Skipping region " << e.region() << " as there are too many haplotypes";
        haplotype_generator.clear_progress();
        return GeneratorStatus::skipped;
    }
    return GeneratorStatus::good;
}

bool Caller::is_saturated(const HaplotypeBlock& haplotypes, const Latents& latents) const
{
    return haplotypes.size() == parameters_.max_haplotypes
           && count_probable_haplotypes(*latents.haplotype_posteriors()) > parameters_.max_haplotypes / 2;
}

unsigned Caller::count_probable_haplotypes(const Caller::Latents::HaplotypeProbabilityMap& haplotype_posteriors) const
{
    return std::count_if(std::cbegin(haplotype_posteriors), std::cend(haplotype_posteriors),
                         [this] (const auto& p) { return p.second >= parameters_.saturation_limit; });
}

void Caller::filter_haplotypes(const bool prefilter_had_removal_impact,
                               const HaplotypeBlock& haplotypes,
                               HaplotypeGenerator& haplotype_generator,
                               const HaplotypeLikelihoodArray& haplotype_likelihoods,
                               const Latents& latents,
                               const std::deque<Haplotype>& protected_haplotypes) const
{
    if (prefilter_had_removal_impact) { // if there was no impact before then there can't be now either
        if (haplotype_generator.removal_has_impact()) {
            const auto max_to_remove = haplotype_generator.max_removal_impact();
            const auto& haplotype_posteriors = *latents.haplotype_posteriors();
            auto removable_haplotypes = get_removable_haplotypes(haplotypes, haplotype_likelihoods,
                                                                 haplotype_posteriors,
                                                                 protected_haplotypes,
                                                                 max_to_remove);
            if (debug_log_) {
                stream(*debug_log_) << "Discarding " << removable_haplotypes.size() << " of "
                                        << max_to_remove << " removable haplotypes with low posterior support";
            }
            haplotype_generator.remove(removable_haplotypes);
        } else if (debug_log_) {
            *debug_log_ << "No posterior haplotype filtering applied as there is no removal impact";
        }
    } else if (debug_log_) {
        *debug_log_ << "No posterior haplotype filtering applied as prefilter had no removal impact";
    }
}

bool Caller::try_early_detect_phase_regions(const MappableBlock<Haplotype>& haplotypes,
                                            const MappableFlatSet<Variant>& candidates,
                                            const GenomicRegion& active_region,
                                            const Latents& latents,
                                            const boost::optional<GenomicRegion>& backtrack_region) const
{
    return parameters_.try_early_phase_detection
        && !backtrack_region
        && count_overlapped(candidates, active_region) > 10
        && haplotypes.size() > parameters_.max_haplotypes / 2;
}

std::vector<decltype(Phaser::PhaseSet::site_indices)>
find_common_phase_regions(Phaser::PhaseSetMap& phasings)
{
    if (phasings.size() > 1) {
        std::map<decltype(Phaser::PhaseSet::site_indices), unsigned> phase_set_sample_assignment_counts {};
        for (auto& p : phasings) {
            for (Phaser::PhaseSet& phase_set : p.second) {
                ++phase_set_sample_assignment_counts[std::move(phase_set.site_indices)];
            }
        }
        const auto num_samples = phasings.size();
        erase_if(phase_set_sample_assignment_counts, [=] (const auto& p) { return p.second < num_samples; });
        return extract_keys(phase_set_sample_assignment_counts);
    } else {
        std::vector<decltype(Phaser::PhaseSet::site_indices)> result {};
        result.reserve(std::cbegin(phasings)->second.size());
        for (auto& phase_set : std::cbegin(phasings)->second) {
            result.push_back(std::move(phase_set.site_indices));
        }
        return result;
    }
}

bool is_contiguous(const decltype(Phaser::PhaseSet::site_indices)& sites)
{
    if (sites.size() < 2) return true;
    const static auto is_discontiguous = [] (auto lhs, auto rhs) { return lhs + 1 < rhs; };
    return std::adjacent_find(std::cbegin(sites), std::cend(sites), is_discontiguous) == std::cend(sites);
}

auto
encompassing_region(const decltype(Phaser::PhaseSet::site_indices)& phase_set,
                    const std::vector<GenomicRegion>& sites)
{
    assert(!phase_set.empty());
    auto result = sites[phase_set.front()];
    std::for_each(std::next(std::cbegin(phase_set)), std::cend(phase_set), [&] (auto site_idx) {
        assert(!begins_before(sites[site_idx], result));
        if (ends_before(result, sites[site_idx])) {
            result = closed_region(result, sites[site_idx]);
        }
    });
    return result;
}

bool 
is_suitable_head_phase_set(const decltype(Phaser::PhaseSet::site_indices)& phase_set,
                           const std::vector<GenomicRegion>& sites)
{
    if (phase_set.front() > 0) return false;
    if (phase_set.back() >= sites.size() - 1) return false;
    if (!is_contiguous(phase_set)) return false;
    const auto phase_region = encompassing_region(phase_set, sites);
    return !(overlaps(phase_region, sites[phase_set.back() + 1]) || are_adjacent(phase_region, sites[phase_set.back() + 1]));
}

boost::optional<GenomicRegion>
Caller::find_phased_head(const MappableBlock<Haplotype>& haplotypes,
                         const MappableFlatSet<Variant>& candidates,
                         const GenomicRegion& active_region,
                         const Latents& latents) const
{
    if (debug_log_) stream(*debug_log_) << "Trying to find complete phase regions in " << active_region;
    const auto active_candidates = contained_range(candidates, active_region);
    const auto viable_phase_regions = extract_regions(active_candidates);
    auto phasings = phaser_.phase(haplotypes, *latents.genotype_posteriors(), viable_phase_regions, get_genotype_calls(latents));
    auto common_phase_regions = find_common_phase_regions(phasings);
    if (common_phase_regions.size() > 1 && is_suitable_head_phase_set(common_phase_regions.front(), viable_phase_regions)) {
        return closed_region(viable_phase_regions.front(), head_region(viable_phase_regions[common_phase_regions.front().back() + 1]));
    } else {
        return boost::none;
    }
}

namespace {

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
                           const HaplotypeBlock& haplotypes,
                           const HaplotypeLikelihoodArray& haplotype_likelihoods,
                           const ReadMap& reads,
                           const Latents& latents,
                           std::deque<CallWrapper>& result,
                           boost::optional<GenomicRegion>& prev_called_region,
                           GenomicRegion& completed_region) const
{
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
        const auto reportable_uncalled_region = overlapped_region(call_region, uncalled_region); // uncalled_region is padded
        if (reportable_uncalled_region) {
            const auto refcall_region = right_overhang_region(*reportable_uncalled_region, completed_region);
            const auto pileups = make_pileups(reads, latents, refcall_region);
            auto alleles = generate_reference_alleles(refcall_region, calls);
            auto reference_calls = call_reference_helper(alleles, latents, pileups);
            const auto itr = utils::append(std::move(reference_calls), calls);
            std::inplace_merge(std::begin(calls), itr, std::end(calls));
        }
    }
    utils::append(std::move(calls), result);
    prev_called_region = uncalled_region;
    completed_region = encompassing_region(completed_region, passed_region);
}

Genotype<IndexedHaplotype<>> Caller::call_genotype(const Latents& latents, const SampleName& sample) const
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
        const auto called_genotype = call_genotype(latents, sample);
        for (const auto& haplotype : called_genotype) {
            result.push_back(haplotype.haplotype());
        }
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
                                  const HaplotypeBlock& haplotypes,
                                  const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    if (parameters_.model_posterior_policy == ModelPosteriorPolicy::all
        || (parameters_.model_posterior_policy == ModelPosteriorPolicy::special && requires_model_evaluation(calls))) {
        const auto mp = calculate_model_posterior(haplotypes, haplotype_likelihoods, latents);
        if (mp) {
            for (auto& call : calls) {
                if (mp->joint) {
                    call->set_model_posterior(probability_false_to_phred(1 - *mp->joint));
                }
                if (!mp->samples.empty()) {
                    for (std::size_t sample_idx {0}; sample_idx < samples_.size(); ++sample_idx) {
                        if (mp->samples[sample_idx]) {
                            call->set_model_posterior(samples_[sample_idx], probability_false_to_phred(1 - *mp->samples[sample_idx]));
                        }
                    }
                }
            }
        }
    }
}

namespace {

void set_phasing(std::vector<CallWrapper>& calls, const Phaser::PhaseSetMap& phasings, const GenomicRegion& calling_region)
{
    for (const auto& p : phasings) {
        const SampleName& sample {p.first};
        for (const Phaser::PhaseSet& phase_set : p.second) {
            const auto phase_set_region = closed_region(calls[phase_set.site_indices.front()],
                                                        calls[phase_set.site_indices.back()]);
            for (const auto call_idx : phase_set.site_indices) {
                assert(call_idx < calls.size());
                calls[call_idx]->set_phase(sample, {phase_set_region, phase_set.quality});
            }
        }
    }
}

} // namespace

namespace debug {

template <typename S, typename Map>
void print_genotype_posteriors(S&& stream,
                               const SampleName& sample,
                               const Map& genotype_posteriors,
                               const std::size_t n)
{
    const auto m = std::min(n, genotype_posteriors.size());
    if (m == genotype_posteriors.size()) {
        stream << "Printing all genotype posteriors for sample " << sample << '\n';
    } else {
        stream << "Printing top " << m << " genotype posteriors for sample " << sample << '\n';
    }
    using GenotypeReference = std::reference_wrapper<const Genotype<IndexedHaplotype<>>>;
    std::vector<std::pair<GenotypeReference, double>> v {};
    v.reserve(genotype_posteriors.size());
    for (const auto& p : genotype_posteriors) {
        v.emplace_back(p.first, p.second);
    }
    const auto mth = std::next(std::begin(v), m);
    std::partial_sort(std::begin(v), mth, std::end(v), [] (const auto& lhs, const auto& rhs) { return lhs.second > rhs.second; });
    std::for_each(std::begin(v), mth, [&] (const auto& p) {
        print_variant_alleles(stream, p.first.get());
        stream << " " << p.second << '\n';
    });
}

template <typename S, typename G>
void print_genotype_posteriors(S&& stream,
                               const ProbabilityMatrix<G>& genotype_posteriors,
                               const std::size_t n = std::numeric_limits<std::size_t>::max())
{
    for (const auto& p : genotype_posteriors) {
        print_genotype_posteriors(stream, p.first, p.second, n);
    }
}

} // namespace debug

void Caller::set_phasing(std::vector<CallWrapper>& calls,
                         const Latents& latents,
                         const HaplotypeBlock& haplotypes,
                         const GenomicRegion& call_region) const
{
    if (debug_log_) stream(*debug_log_) << "Phasing " << calls.size() << " calls in " << call_region;
    if (trace_log_) debug::print_genotype_posteriors(stream(*trace_log_), *latents.genotype_posteriors());
    const auto call_regions = extract_regions(calls);
    const auto phase_sets = phaser_.phase(haplotypes, *latents.genotype_posteriors(), call_regions, get_genotype_calls(latents));
    if (debug_log_) debug::print_phase_sets(stream(*debug_log_), phase_sets, call_regions);
    octopus::set_phasing(calls, phase_sets, call_region);
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

MappableFlatSet<Variant> 
Caller::generate_candidate_variants(const GenomicRegion& region, OptionalThreadPool workers) const
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

HaplotypeGenerator 
Caller::make_haplotype_generator(const MappableFlatSet<Variant>& candidates,
                                 const ReadMap& reads, 
                                 const boost::optional<TemplateMap>& read_templates) const
{
    if (read_templates) {
        return haplotype_generator_builder_.build(reference_, candidates, reads, *read_templates);
    } else {
        return haplotype_generator_builder_.build(reference_, candidates, reads, boost::none);
    }
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

HaplotypeLikelihoodArray::FlankState
calculate_flank_state(const MappableBlock<Haplotype>& haplotypes,
                      const GenomicRegion& active_region,
                      const MappableFlatSet<Variant>& candidates)
{
    const auto flank_regions = calculate_flank_regions(mapped_region(haplotypes), active_region, candidates);
    return {size(flank_regions.first), size(flank_regions.second)};
}

bool Caller::compute_haplotype_likelihoods(HaplotypeLikelihoodArray& haplotype_likelihoods,
                                           const GenomicRegion& active_region,
                                           const HaplotypeBlock& haplotypes,
                                           const MappableFlatSet<Variant>& candidates,
                                           const boost::variant<ReadMap, TemplateMap>& active_reads,
                                           OptionalThreadPool workers) const
{
    assert(haplotype_likelihoods.is_empty());
    boost::optional<HaplotypeLikelihoodArray::FlankState> flank_state {};
    if (debug_log_) {
        stream(*debug_log_) << "Calculating likelihoods for " << haplotypes.size() << " haplotypes";
        debug::print_active_candidates(stream(*debug_log_), candidates, active_region);
    }
    if (likelihood_model_.can_use_flank_state()) {
        flank_state = calculate_flank_state(haplotypes, active_region, candidates);
        if (debug_log_) {
            debug::print_inactive_flanking_candidates(stream(*debug_log_), candidates, active_region,
                                                      mapped_region(haplotypes));
        }
    }
    try {
        boost::apply_visitor([&] (const auto& reads) { 
            haplotype_likelihoods.populate(reads, haplotypes, std::move(flank_state), workers); }, active_reads);
    } catch(const HaplotypeLikelihoodModel::ShortHaplotypeError& e) {
        if (debug_log_) {
            stream(*debug_log_) << "Skipping " << active_region << " as a haplotype was too short by "
                                << e.required_extension() << "bp";
        }
        return false;
    }
    if (trace_log_) {
        boost::apply_visitor([&] (const auto& reads) { 
            debug::print_read_haplotype_likelihoods(stream(*trace_log_), haplotypes, reads,
                                               haplotype_likelihoods); }, active_reads);
    }
    return true;
}

std::vector<Haplotype>
Caller::filter(HaplotypeBlock& haplotypes,
               const HaplotypeLikelihoodArray& haplotype_likelihoods,
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
Caller::get_removable_haplotypes(const HaplotypeBlock& haplotypes,
                                 const HaplotypeLikelihoodArray& haplotype_likelihoods,
                                 const Caller::Latents::HaplotypeProbabilityMap& haplotype_posteriors,
                                 const std::deque<Haplotype>& protected_haplotypes,
                                 const unsigned max_to_remove) const
{
    if (debug_log_) {
        stream(*debug_log_) << "Protecting " << protected_haplotypes.size() << " haplotypes from filtering";
    }
    HaplotypeReferenceProbabilityMap haplotype_ref_posteriors {};
    haplotype_ref_posteriors.reserve(haplotype_posteriors.size());
    for (const auto& p : haplotype_posteriors) {
        haplotype_ref_posteriors.emplace(p.first.haplotype(), p.second);
    }
    std::vector<std::reference_wrapper<const Haplotype>> result {};
    if (protected_haplotypes.empty()) {
        result = extract_removable(haplotypes, haplotype_ref_posteriors, samples_, haplotype_likelihoods,
                                 max_to_remove, parameters_.haplotype_extension_threshold);
    } else if (protected_haplotypes.size() == 1
            && parameters_.protect_reference_haplotype
            && is_reference(protected_haplotypes.front())) {
        result = extract_removable(haplotypes, haplotype_ref_posteriors, samples_, haplotype_likelihoods,
                                        max_to_remove, parameters_.haplotype_extension_threshold);
        auto reference_itr = find_reference(result);
        if (reference_itr != std::cend(result)) result.erase(reference_itr);
    } else {
        std::vector<Haplotype> removable_haplotypes {};
        removable_haplotypes.reserve(haplotypes.size());
        assert(std::is_sorted(std::cbegin(haplotypes), std::cend(haplotypes)));
        assert(std::is_sorted(std::cbegin(protected_haplotypes), std::cend(protected_haplotypes)));
        std::set_difference(std::cbegin(haplotypes), std::cend(haplotypes),
                            std::cbegin(protected_haplotypes), std::cend(protected_haplotypes),
                            std::back_inserter(removable_haplotypes));
        result = extract_removable(removable_haplotypes, haplotype_ref_posteriors, samples_, haplotype_likelihoods,
                                 max_to_remove, parameters_.haplotype_extension_threshold);
    }
    return result;
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

auto get_safe_haplotype_region(const GenomicRegion& region, const HaplotypeLikelihoodModel& model, const ReadMap& reads)
{
    if (has_coverage(reads)) {
        auto min_safe_pad = model.pad_requirement() + (max_read_length(reads) / 2) + 1;
        return expand(encompassing_region(reads), min_safe_pad);
    } else {
        return region;
    }
}

std::vector<CallWrapper> Caller::call_reference(const GenomicRegion& region, const ReadMap& reads) const
{
    const auto active_reads = copy_overlapped(reads, region);
    const HaplotypeBlock haplotypes {Haplotype {get_safe_haplotype_region(region, likelihood_model_, active_reads), reference_}};
    auto haplotype_likelihoods = make_haplotype_likelihood_cache();
    haplotype_likelihoods.populate(active_reads, haplotypes);
    const auto latents = infer_latents(haplotypes, haplotype_likelihoods);
    const auto pileups = make_pileups(active_reads, *latents, region);
    const auto alleles = generate_reference_alleles(region);
    return call_reference_helper(alleles, *latents, pileups);
}

std::vector<CallWrapper>
Caller::call_reference_helper(const std::vector<Allele>& alleles, const Latents& latents, const ReadPileupMap& pileups) const
{
    auto refcalls = call_reference(alleles, latents, pileups);
    if (parameters_.max_refcall_posterior) {
        for (auto& refcall : refcalls) {
            refcall->set_quality(std::min(refcall->quality(), *parameters_.max_refcall_posterior));
        }
    }
    if (is_merge_block_refcalling()) {
        refcalls = squash_reference_calls(std::move(refcalls));
    }
    return wrap(std::move(refcalls));
}

namespace {

template <typename Container>
auto extract_unique_regions(const Container& mappables)
{
    auto result = extract_regions(mappables);
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

auto extract_called_regions(const std::vector<CallWrapper>& calls)
{
    auto regions = extract_regions(calls);
    const auto is_insertion = [] (const auto& region) { return is_empty(region); };
    auto insertion_itr = std::find_if(std::begin(regions), std::end(regions), is_insertion);
    if (insertion_itr != std::end(regions)) {
        do {
            *insertion_itr = expand_lhs(*insertion_itr, 1);
            insertion_itr = std::find_if(std::next(insertion_itr), std::end(regions), is_insertion);
        } while (insertion_itr != std::end(regions));
        std::sort(std::begin(regions), std::end(regions));
    }
    return extract_covered_regions(regions);
}

auto extract_uncalled_reference_regions(const GenomicRegion& region,
                                        const std::vector<CallWrapper>& calls)
{
    auto called_regions = extract_called_regions(calls);
    return extract_intervening_regions(overlap_range(called_regions, region), region);
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
                                   const std::vector<CallWrapper>& calls) const
{
    auto refcall_regions = extract_uncalled_reference_regions(region, calls);
    if (parameters_.refcall_type == RefCallType::positional || is_merge_block_refcalling()) {
        return make_positional_reference_alleles(std::move(refcall_regions), reference_);
    } else {
        return make_reference_alleles(std::move(refcall_regions), reference_);
    }
}

std::vector<Allele> Caller::generate_reference_alleles(const GenomicRegion& region) const
{
    return generate_reference_alleles(region, {});
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
        const auto called_genotype = genotype_cast<Haplotype>(call_genotype(latents, sample));
        result.emplace(sample, octopus::make_pileups(reads.at(sample), called_genotype, region));
    }
    return result;
}

namespace {

std::unique_ptr<ReferenceCall>
squash(const std::vector<std::unique_ptr<ReferenceCall>>& refcalls, 
       const std::vector<SampleName>& samples,
       const Phred<double> quality)
{
    assert(!refcalls.empty());
    auto region = encompassing_region(refcalls.front()->mapped_region(), refcalls.back()->mapped_region());
    Allele::NucleotideSequence sequence {};
    sequence.reserve(size(region));
    for (const auto& refcall : refcalls) {
        sequence.insert(std::cend(sequence), std::cbegin(refcall->reference().sequence()), std::cend(refcall->reference().sequence()));
    }
    Allele reference {std::move(region), std::move(sequence)};
    std::map<SampleName, ReferenceCall::GenotypeCall> genotypes {};
    for (const auto& sample : samples) {
        const auto& genotype_call = refcalls.front()->get_genotype_call(sample);
        genotypes.emplace(sample, ReferenceCall::GenotypeCall {genotype_call.genotype.ploidy(), genotype_call.posterior});
    }
    return std::make_unique<ReferenceCall>(std::move(reference), quality, std::move(genotypes));
}

} // namespace

std::vector<std::unique_ptr<ReferenceCall>>
Caller::squash_reference_calls(std::vector<std::unique_ptr<ReferenceCall>> refcalls) const
{
    assert(parameters_.refcall_block_merge_threshold);
    if (refcalls.size() < 2) return refcalls;
    // Implements a simple greedy hierarchal clustering algorithm, with the
    // restriction only adjacent clusters (i.e. contiguous runs of refcalls) can be joined.
    using ClusterID = std::size_t;
    std::vector<ClusterID> position_cluster_assignments(refcalls.size());
    // Initialise assignments
    std::iota(std::begin(position_cluster_assignments), std::end(position_cluster_assignments), ClusterID {0});
    using IndexRangePair = std::pair<std::size_t, std::size_t>; 
    std::unordered_map<ClusterID, std::pair<IndexRangePair, double>> clusters {};
    clusters.reserve(refcalls.size() + 1);
    // Initialise clusters
    for (ClusterID id {0}; id < refcalls.size(); ++id) {
        clusters.emplace(id, std::make_pair(std::make_pair(id, id), refcalls[id]->quality().score()));
    }
    std::priority_queue<std::pair<double, std::pair<ClusterID, ClusterID>>> candidates {};
    // Initialise candidates
    for (ClusterID lhs_id {0}; lhs_id < (refcalls.size() - 1); ++lhs_id) {
        if (are_adjacent(*refcalls[lhs_id], *refcalls[lhs_id + 1])) {
            const auto qual_diff = std::abs(refcalls[lhs_id]->quality().score() - refcalls[lhs_id + 1]->quality().score());
            if (qual_diff <= parameters_.refcall_block_merge_threshold->score()) {
                candidates.emplace(qual_diff, std::make_pair(lhs_id, lhs_id + 1));
            }
        }
    }
    std::size_t next_cluster_id = refcalls.size();
    while (!candidates.empty()) {
        const auto candidate = candidates.top();
        candidates.pop();
        const auto lhs_cluster_itr = clusters.find(candidate.second.first);
        const auto rhs_cluster_itr = clusters.find(candidate.second.second);
        // If we don't find a cluster then it must have already been merged into another cluster
        if (lhs_cluster_itr != std::cend(clusters) && rhs_cluster_itr != std::cend(clusters)) {
            auto new_cluster_range = std::make_pair(lhs_cluster_itr->second.first.first, rhs_cluster_itr->second.first.second);
            std::vector<double> qualities(new_cluster_range.second - new_cluster_range.first + 1);
            std::transform(std::next(std::cbegin(refcalls), new_cluster_range.first),
                           std::next(std::cbegin(refcalls), new_cluster_range.second + 1),
                           std::begin(qualities),
                           [] (const auto& call) { return call->quality().score(); });
            const auto new_cluster_quality = maths::median(qualities);
            clusters.emplace(next_cluster_id, std::make_pair(new_cluster_range, new_cluster_quality));
            std::fill(std::next(std::begin(position_cluster_assignments), new_cluster_range.first),
                      std::next(std::begin(position_cluster_assignments), new_cluster_range.second + 1),
                      next_cluster_id);
            if (new_cluster_range.first > 0
             && are_adjacent(*refcalls[new_cluster_range.first - 1], *refcalls[new_cluster_range.first])) {
                const auto lhs_cluster_id = position_cluster_assignments[new_cluster_range.first - 1];
                const auto lhs_cluster_qual = clusters.at(lhs_cluster_id).second;
                const auto lhs_qual_diff = std::abs(new_cluster_quality - lhs_cluster_qual);
                if (lhs_qual_diff <= parameters_.refcall_block_merge_threshold->score()) {
                    candidates.emplace(lhs_qual_diff, std::make_pair(lhs_cluster_id, next_cluster_id));
                }
            }
            if (new_cluster_range.second < (refcalls.size() - 1)
             && are_adjacent(*refcalls[new_cluster_range.second], *refcalls[new_cluster_range.second + 1])) {
                const auto rhs_cluster_id = position_cluster_assignments[new_cluster_range.second + 1];
                const auto rhs_cluster_qual = clusters.at(rhs_cluster_id).second;
                const auto rhs_qual_diff = std::abs(new_cluster_quality - rhs_cluster_qual);
                if (rhs_qual_diff <= parameters_.refcall_block_merge_threshold->score()) {
                    candidates.emplace(rhs_qual_diff, std::make_pair(next_cluster_id, rhs_cluster_id));
                }
            }
            ++next_cluster_id;
            clusters.erase(lhs_cluster_itr);
            clusters.erase(rhs_cluster_itr);
        }
    }
    std::vector<std::unique_ptr<ReferenceCall>> result {};
    result.reserve(clusters.size());
    for (const auto& p : clusters) {
        std::size_t lhs_idx {p.second.first.first}, rhs_idx {p.second.first.second};
        std::vector<std::unique_ptr<ReferenceCall>> cluster(rhs_idx - lhs_idx + 1);
        std::move(std::next(std::begin(refcalls), lhs_idx),
                  std::next(std::begin(refcalls), rhs_idx + 1),
                  std::begin(cluster));
        result.push_back(squash(cluster, samples_, Phred<> {p.second.second}));
    }
    std::sort(std::begin(result), std::end(result), [] (const auto& lhs, const auto& rhs) { return *lhs < *rhs; });
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

template <typename S, typename Container>
void print_haplotypes(S&& stream, const Container& haplotypes,
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

template <typename Container>
void print_haplotypes(const Container& haplotypes, const Resolution resolution)
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
    for (const auto& p : haplotype_posteriors) {
        v.emplace_back(p.first.haplotype(), p.second);
    }
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
