// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "indel_profiler.hpp"

#include <iterator>
#include <algorithm>
#include <utility>
#include <thread>
#include <cassert>
#include <iostream>

#include "basics/genomic_region.hpp"
#include "basics/cigar_string.hpp"
#include "utils/genotype_reader.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_header.hpp"
#include "utils/append.hpp"
#include "utils/repeat_finder.hpp"
#include "utils/read_stats.hpp"
#include "read_realigner.hpp"

#include "logging/progress_meter.hpp"

namespace octopus {

namespace {

unsigned get_pool_size(const IndelProfiler::PerformanceConfig& config)
{
    const auto num_cores = std::thread::hardware_concurrency();
    if (config.max_threads) {
        if (*config.max_threads > 1) {
            return num_cores > 0 ? std::min(*config.max_threads, num_cores) : *config.max_threads;
        } else {
            return 0;
        }
    } else {
        return num_cores > 0 ? num_cores : 8;
    }
}

} // namespace

IndelProfiler::IndelProfiler(ProfileConfig config)
: config_ {std::move(config)}
, performance_config_ {}
, workers_ {get_pool_size(performance_config_)}
{}

IndelProfiler::IndelProfiler(ProfileConfig config, PerformanceConfig performance_config)
: config_ {std::move(config)}
, performance_config_ {}
, workers_ {get_pool_size(performance_config_)}
{}

IndelProfiler::IndelProfile
IndelProfiler::profile(const ReadPipe& src, VcfReader& variants, const ReferenceGenome& reference) const
{
    const auto samples = src.samples();
    check_samples(samples, variants);
    IndelProfile result {};
    boost::optional<GenomicRegion> batch_region {};
    GenomicRegion analysis_region;
    for (auto p = variants.iterate(); !batch_region || p.first != p.second; ) {
        if (!batch_region || !is_same_contig(analysis_region, *p.first)) {
            analysis_region = reference.contig_region(contig_name(*p.first));
        }
        const auto data = read_next_data_batch(p.first, p.second, src, reference, samples, analysis_region, batch_region);
        evaluate_indel_profile(data, result);
        batch_region = mapped_region(data.reference);
    }
    return result;
}

IndelProfiler::IndelProfile
IndelProfiler::profile(const ReadPipe& src, VcfReader& variants, const ReferenceGenome& reference,
                       const InputRegionMap& regions) const
{
    const auto samples = src.samples();
    check_samples(samples, variants);
    ProgressMeter progress {regions};
    progress.start();
    IndelProfile result {};
    for (const auto& r : regions) {
        for (const auto& analysis_region : r.second) {
            boost::optional<GenomicRegion> batch_region {};
            for (auto p = variants.iterate(analysis_region); !batch_region || p.first != p.second; ) {
                const auto data = read_next_data_batch(p.first, p.second, src, reference, samples, analysis_region, batch_region);
                evaluate_indel_profile(data, result);
                batch_region = mapped_region(data.reference);
                progress.log_completed(*batch_region);
            }
        }
    }
    progress.stop();
    return result;
}

IndelProfiler::IndelProfile
IndelProfiler::profile(const ReadPipe& src, VcfReader& variants, const ReferenceGenome& reference,
                       const GenomicRegion& region) const
{
    const auto samples = src.samples();
    check_samples(samples, variants);
    IndelProfile result {};
    boost::optional<GenomicRegion> batch_region {};
    for (auto p = variants.iterate(region); !batch_region || p.first != p.second; ) {
        const auto data = read_next_data_batch(p.first, p.second, src, reference, samples, region, batch_region);
        evaluate_indel_profile(data, result);
        batch_region = mapped_region(data.reference);
    }
    return result;
}

// private methods

void IndelProfiler::check_samples(const SampleList& samples, const VcfReader& variants) const
{
    auto vcf_samples = variants.fetch_header().samples();
    std::sort(std::begin(vcf_samples), std::end(vcf_samples));
    if (!std::all_of(std::cbegin(samples), std::cend(samples), [&] (const auto& sample) {
                     return std::binary_search(std::cbegin(vcf_samples), std::cend(vcf_samples), sample); })) {
        throw std::runtime_error {"Not all samples in VCF"};
    }
}

namespace {

GenomicRegion get_phase_set(const VcfRecord& record, const SampleName& sample)
{
    auto result = get_phase_region(record, sample);
    return result ? *result : mapped_region(record);
}

std::vector<GenomicRegion> get_phase_sets(const VcfRecord& record, const std::vector<SampleName>& samples)
{
    std::vector<GenomicRegion> result{};
    result.reserve(samples.size());
    std::transform(std::cbegin(samples), std::cend(samples), std::back_inserter(result),
                   [&record](const auto& sample) { return get_phase_set(record, sample); });
    return result;
}

GenomicRegion get_phase_region(const VcfRecord& record, const std::vector<SampleName>& samples)
{
    return encompassing_region(get_phase_sets(record, samples));
}

template<typename T, typename _>
std::vector<T> copy_each_first(const std::vector<std::pair<T, _>>& items)
{
    std::vector<T> result {};
    result.reserve(items.size());
    std::transform(std::cbegin(items), std::cend(items), std::back_inserter(result),
                   [] (const auto& p) { return p.first; });
    return result;
}

} // namespace

IndelProfiler::CallBlock
IndelProfiler::read_next_block(VcfIterator& first, const VcfIterator& last, const SampleList& samples) const
{
    std::vector<std::pair<VcfRecord, GenomicRegion>> block {};
    for (; first != last; ++first) {
        const VcfRecord& call {*first};
        auto call_phase_region = get_phase_region(call, samples);
        if (!block.empty() && !overlaps(block.back().second, call_phase_region)
            && inner_distance(block.back().second, call_phase_region) > 2 * config_.max_length) {
            return copy_each_first(block);
        }
        block.emplace_back(call, std::move(call_phase_region));
    }
    return copy_each_first(block);
}

namespace {

Genotype<Haplotype> remap(const Genotype<Haplotype>& genotype, const GenomicRegion& region)
{
    Genotype<Haplotype> result {genotype.ploidy()};
    for (const auto& haplotype : genotype) {
        result.emplace(remap(haplotype, region));
    }
    return result;
}

} // namespace

IndelProfiler::DataBatch
IndelProfiler::read_next_data_batch(VcfIterator& first, const VcfIterator& last, const ReadPipe& src,
                                    const ReferenceGenome& reference, const SampleList& samples,
                                    const GenomicRegion& analysis_region,
                                    const boost::optional<GenomicRegion>& prev_batch_region) const
{
    const auto records = read_next_block(first, last, samples);
    GenomicRegion batch_region;
    if (!records.empty()) {
        batch_region = encompassing_region(records);
        if (prev_batch_region && is_same_contig(batch_region, *prev_batch_region)) {
            batch_region = expand_lhs(batch_region, intervening_region_size(*prev_batch_region, batch_region) / 2);
        } else {
            batch_region = closed_region(analysis_region, batch_region);
        }
        if (first != last && is_same_contig(batch_region, *first) && overlaps(*first, analysis_region)) {
            batch_region = expand_rhs(batch_region, intervening_region_size(batch_region, *first) / 2);
        } else {
            batch_region = closed_region(batch_region, analysis_region);
        }
    } else if (prev_batch_region) {
        batch_region = right_overhang_region(analysis_region, *prev_batch_region);
    } else {
        batch_region = analysis_region;
    }
    DataBatch result {Haplotype {std::move(batch_region), reference}, {}, {}, {}};
    if (is_empty_region(result.reference)) return result;
    auto reads = src.fetch_reads(mapped_region(result.reference));
    if (has_coverage(reads)) {
        result.repeats[result.reference] = find_repeats(result.reference);
    } else {
        const auto reads_region = encompassing_region(reads);
        if (contains(result.reference, reads_region)) {
            result.repeats[result.reference] = find_repeats(result.reference);
        } else {
            const auto expanded_reference = remap(result.reference, reads_region);
            result.repeats[result.reference] = find_repeats(expanded_reference);
        }
    }
    if (!records.empty()) {
        auto genotypes = extract_genotypes(records, samples, reference);
        result.support.reserve(2 * samples.size()); // guess
        for (const auto& sample : samples) {
            if (prev_batch_region) reads[sample].erase_overlapped(*prev_batch_region);
            auto& sample_genotypes = genotypes.at(sample);
            for (auto& genotype : sample_genotypes) {
                genotype = remap(genotype, expand(mapped_region(genotype), 2 * config_.max_length));
            }
            Genotype<Haplotype> sample_genotype;
            if (sample_genotypes.size() == 1) {
                sample_genotype = std::move(sample_genotypes.front());
            } else {
                for (auto& genotype : sample_genotypes) {
                    for (auto& haplotype : genotype) sample_genotype.emplace(haplotype);
                }
            }
            evaluate_support(sample_genotype, reads[sample], result.support);
            reads[sample].clear();
            reads[sample].shrink_to_fit();
            for (const auto& record : records) {
                std::vector<Allele> alleles; bool has_ref;
                std::tie(alleles, has_ref) = get_called_alleles(record, sample);
                if (has_ref) alleles.erase(std::cbegin(alleles));
                alleles.erase(std::remove_if(std::begin(alleles), std::end(alleles),
                                             [] (const auto& allele) { return !is_indel(allele); }),
                              std::end(alleles));
                result.indels.insert(std::begin(alleles), std::end(alleles));
            }
        }
    } else {
        const Genotype<Haplotype> reference_genotype {result.reference};
        for (auto& p : reads) {
            evaluate_support(reference_genotype, reads[p.first], result.support);
            p.second.clear();
            p.second.shrink_to_fit();
        }
    }
    result.repeats.reserve(result.support.size() + 1);
    for (auto& p : result.support) {
        if (result.repeats.count(p.first) == 0) {
            result.repeats[p.first] = find_repeats(p.first);
        }
        p.second.shrink_to_fit();
    }
    return result;
}

namespace {

boost::optional<Allele> find_indel_error(const AlignedRead& read, const GenomicRegion& region)
{
    if (!contains(read, region) || !has_indel(read)) return boost::none;
    auto read_section = copy_sequence(read, mapped_region(region));
    if (read_section.size() != size(region)) {
        auto error_region = region;
        if (read_section.size() < size(region)) {
            error_region = expand_rhs(error_region, -static_cast<GenomicRegion::Distance>(read_section.size()));
            read_section.clear();
        } else {
            read_section.erase(std::prev(std::cend(read_section), size(region)), std::cend(read_section));
            error_region = head_region(error_region);
        }
        return Allele {std::move(error_region), std::move(read_section)};
    } else {
        return boost::none;
    }
}

bool is_closest_repeat(const TandemRepeat& repeat, const Allele& indel, const MappableFlatSet<TandemRepeat>& repeats)
{
    const auto candidates = overlap_range(repeats, indel);
    if (size(candidates) > 1) {
        const auto matches_indel_begin = [&] (const auto& repeat) { return begins_equal(repeat, indel); };
        const auto perfect_candidate_itr = std::find_if(std::cbegin(candidates), std::cend(candidates), matches_indel_begin);
        if (perfect_candidate_itr != std::cend(candidates)) {
            return repeat == *perfect_candidate_itr;
        } else {
            return repeat == *max_overlapped(repeats, indel);
        }
    }
    return true;
}

IndelProfiler::IndelProfile::RepeatState&
find_or_insert_motif(const Haplotype::NucleotideSequence& motif, std::deque<IndelProfiler::IndelProfile::RepeatState>& states)
{
    const auto motif_less = [] (const auto& lhs, const auto& rhs) { return lhs.motif < rhs; };
    auto motif_itr = std::lower_bound(std::begin(states), std::end(states), motif, motif_less);
    if (motif_itr != std::end(states)) {
        if (motif_itr->motif == motif) {
            return *motif_itr;
        }
    }
    IndelProfiler::IndelProfile::RepeatState state {};
    state.motif = motif;
    return *states.insert(motif_itr, std::move(state));
}

} // namespace

void IndelProfiler::evaluate_indel_profile(const DataBatch& data, IndelProfile& result) const
{
    if (is_empty_region(data.reference)) return;
    assert(data.repeats.count(data.reference) == 1);
    const auto& reference_repeats = data.repeats.at(data.reference);
    // Add polymorphisms
    for (const auto& repeat : reference_repeats) {
        if (result.states.size() <= repeat.period()) result.states.resize(repeat.period() + 1);
        const auto periods = count_periods(repeat);
        if (result.states[repeat.period()].size() <= periods) result.states[repeat.period()].resize(periods + 1);
        IndelProfile::RepeatState& state {find_or_insert_motif(repeat.motif(), result.states[repeat.period()][periods])};
        state.span += region_size(repeat);
        ++state.reference_count;
        for (const auto& indel : overlap_range(data.indels, repeat)) {
            if (is_closest_repeat(repeat, indel, reference_repeats)) {
                const auto indel_size = reference_distance(indel);
                if (state.polymorphism_counts.size() <= indel_size) state.polymorphism_counts.resize(indel_size + 1);
                ++state.polymorphism_counts[indel_size];
            }
        }
    }
    // Add errors
    for (const auto& p : data.support) {
        const Haplotype& haplotype {p.first};
        const auto& haplotype_repeats = data.repeats.at(haplotype);
        for (const auto& read : p.second) {
            for (const auto& repeat : contained_range(haplotype_repeats, read)) {
                if (result.states.size() <= repeat.period()) result.states.resize(repeat.period() + 1);
                const auto periods = count_periods(repeat);
                if (result.states[repeat.period()].size() <= periods) result.states[repeat.period()].resize(periods + 1);
                IndelProfile::RepeatState& state {find_or_insert_motif(repeat.motif(), result.states[repeat.period()][periods])};
                ++state.read_count;
                const auto indel_error = find_indel_error(read, mapped_region(repeat));
                if (indel_error && is_closest_repeat(repeat, *indel_error, haplotype_repeats)) {
                    const auto indel_error_length = reference_distance(*indel_error);
                    if (state.error_counts.size() <= indel_error_length) state.error_counts.resize(indel_error_length + 1);
                    ++state.error_counts[indel_error_length];
                }
            }
        }
    }
    const auto non_repeat_regions = extract_intervening_regions(extract_covered_regions(reference_repeats), data.reference);
    for (const auto& region : non_repeat_regions) {
        if (result.states.empty()) result.states.resize(1);
        if (result.states[0].empty()) result.states[0].resize(1);
        if (result.states[0][0].empty()) result.states[0][0].resize(1);
        IndelProfile::RepeatState& state {result.states[0][0][0]};
        state.motif = config_.complex_motif;
        state.reference_count = 1;
        state.span += size(region);
        // Add polymorphisms
        for (const auto& indel : contained_range(data.indels, region)) {
            if (!has_overlapped(data.repeats, indel)) {
                const auto indel_size = reference_distance(indel);
                if (state.polymorphism_counts.size() <= indel_size) state.polymorphism_counts.resize(indel_size + 1);
                ++state.polymorphism_counts[indel_size];
            }
        }
        // Add errors
        for (const auto& p : data.support) {
            for (const auto& read : overlap_range(p.second, region)) {
                ++state.read_count;
                const auto indel_error = find_indel_error(read, region);
                if (indel_error) {
                    const auto indel_error_length = reference_distance(*indel_error);
                    if (state.error_counts.size() <= indel_error_length) state.error_counts.resize(indel_error_length + 1);
                    ++state.error_counts[indel_error_length];
                }
            }
        }
    }
}

namespace {

double error_expectation(const AlignedRead::BaseQualityVector& qualities)
{
    return std::accumulate(std::cbegin(qualities), std::cend(qualities), 0.0,
                           [] (auto curr, auto q) { return curr + std::pow(10.0, -q / 10); });
}

double error_expectation(const AlignedRead& read)
{
    return error_expectation(read.base_qualities());
}

unsigned count_errors(const CigarString& cigar)
{
    return std::accumulate(std::cbegin(cigar), std::cend(cigar), 0u,
                           [] (auto curr, auto op) { return curr + (is_match(op) ? 0u : op.size()); });
}

unsigned count_errors(const AlignedRead& read)
{
    return count_errors(read.cigar());
}

bool over_hmm_band_limit(const CigarString& cigar) noexcept
{
    auto gap_len = static_cast<unsigned>(std::max(std::abs(max_indel_size(cigar)), std::abs(sum_indel_sizes(cigar))));
    return gap_len > HaplotypeLikelihoodModel::pad_requirement();
}

bool is_likely_misaligned(const AlignedRead& read)
{
    if (over_hmm_band_limit(read.cigar())) return true;
    const auto expected_errors = static_cast<unsigned>(std::ceil(error_expectation(read)));
    const auto max_expected_errors = 5 * expected_errors;
    const auto observed_errors = count_errors(read);
    return observed_errors > max_expected_errors;
}

} // namespace

void IndelProfiler::evaluate_support(const Genotype<Haplotype>& genotype, ReadContainer& reads,
                                     HaplotypeReadSupportMap& result) const
{
    // TODO: Need positional gap extension penalties in read likelihood calculation.
    assert(genotype.ploidy() > 0);
    std::vector<AlignedRead> buffer {std::make_move_iterator(std::begin(reads)), std::make_move_iterator(std::end(reads))};
    HaplotypeSupportMap haplotype_support;
    if (genotype.ploidy() > 1) {
        const AssignmentConfig assigner_config {AssignmentConfig::AmbiguousAction::random};
        haplotype_support = compute_haplotype_support(genotype, buffer, assigner_config);
    } else {
        haplotype_support.emplace(genotype[0], std::move(buffer));
    }
    buffer.clear();
    buffer.shrink_to_fit();
    for (auto& p : haplotype_support) {
        safe_realign(p.second, p.first);
        if (config_.check_read_misalignments) {
            p.second.erase(std::remove_if(std::begin(p.second), std::end(p.second),
                                          [] (const auto& read) { return is_likely_misaligned(read); }),
                           std::end(p.second));
            p.second.shrink_to_fit();
        }
        std::sort(std::begin(p.second), std::end(p.second));
        result[p.first].insert(std::make_move_iterator(std::begin(p.second)), std::make_move_iterator(std::end(p.second)));
    }
}

namespace {

void remap_to_reference(TandemRepeat& repeat, const GenomicRegion& origin_region, const CigarString& origin_reference_cigar)
{
    const auto repeat_offset = begin_distance(origin_region, repeat);
    assert(repeat_offset >= 0);
    const auto cigar_before_repeat = copy_sequence(origin_reference_cigar, 0, repeat_offset);
    const auto indel_length = sum_indel_sizes(cigar_before_repeat);
    repeat.region() = shift(repeat.region(), -indel_length);
}

void remap_to_reference(std::vector<TandemRepeat>& repeats, const Haplotype& origin)
{
    const auto& origin_region = mapped_region(origin);
    const auto origin_reference_cigar = origin.cigar();
    if (!has_indel(origin_reference_cigar)) return;
    for (auto& repeat : repeats) remap_to_reference(repeat, origin_region, origin_reference_cigar);
}

auto find_tandem_repeats(const Haplotype& haplotype, const unsigned min_period, unsigned max_period)
{
    auto result = find_exact_tandem_repeats(haplotype.sequence(), mapped_region(haplotype), min_period, max_period);
    if (!is_reference(haplotype)) {
        remap_to_reference(result, haplotype);
    }
    return result;
}

bool is_maximal(const TandemRepeat& repeat, const MappableFlatSet<TandemRepeat>& repeats)
{
    const auto overlapped = overlap_range(repeats, repeat);
    const auto is_contained_by = [&] (const auto& other) { return contains(other, repeat); };
    return std::count_if(std::cbegin(overlapped), std::cend(overlapped), is_contained_by) == 1;
}

MappableFlatSet<TandemRepeat> find_minimum_spanning_set(const MappableFlatSet<TandemRepeat>& repeats)
{
    MappableFlatSet<TandemRepeat> result {};
    std::copy_if(std::cbegin(repeats), std::cend(repeats), std::inserter(result, std::begin(result)),
                 [&] (const auto& repeat) { return is_maximal(repeat, repeats); });
    return result;
}

MappableFlatSet<TandemRepeat> find_minimum_spanning_set(std::vector<TandemRepeat> repeats)
{
    const MappableFlatSet<TandemRepeat> unique_repeats {std::make_move_iterator(std::begin(repeats)),
                                                        std::make_move_iterator(std::end(repeats))};
    return find_minimum_spanning_set(unique_repeats);
}

} // namespace

MappableFlatSet<TandemRepeat> IndelProfiler::find_repeats(const Haplotype& haplotype) const
{
    auto repeats = find_tandem_repeats(haplotype, config_.min_period, config_.max_period);
    repeats.erase(std::remove_if(std::begin(repeats), std::end(repeats),
                                 [this] (const auto& repeat) {
                                     const auto periods = count_periods(repeat);
                                     return periods < config_.min_periods || periods > config_.max_periods
                                          || region_size(repeat) > config_.max_length;
                                 }), std::end(repeats));
    return find_minimum_spanning_set(std::move(repeats));
}

IndelProfiler::IndelProfile
profile_indels(const ReadPipe& reads, VcfReader::Path variants, const ReferenceGenome& reference)
{
    VcfReader vcf {std::move(variants)};
    IndelProfiler profiler {};
    return profiler.profile(reads, vcf, reference);
}

IndelProfiler::IndelProfile
profile_indels(const ReadPipe& reads, VcfReader::Path variants, const ReferenceGenome& reference, const InputRegionMap& regions)
{
    VcfReader vcf {std::move(variants)};
    IndelProfiler profiler {};
    return profiler.profile(reads, vcf, reference, regions);
}

IndelProfiler::IndelProfile
profile_indels(const ReadPipe& reads, VcfReader::Path variants, const ReferenceGenome& reference, const GenomicRegion& region)
{
    VcfReader vcf {std::move(variants)};
    IndelProfiler profiler {};
    return profiler.profile(reads, vcf, reference, region);
}

std::ostream& operator<<(std::ostream& os, const IndelProfiler::IndelProfile::RepeaStateArray& states)
{
    os << "period,periods,motif,reference_count,reference_span,indel_length,polymorphisms,errors,reads";
    if (states.empty()) return os;
    const auto max_period = states.size() - 1;
    for (std::size_t period {0}; period <= max_period; ++period) {
        if (states[period].empty()) continue;
        const auto max_periods = states[period].size() - 1;
        for (std::size_t periods {0}; periods <= max_periods; ++periods) {
            if (!((period == 0 && periods == 0) || (period >= 1 && periods > 1))) continue;
            using State = IndelProfiler::IndelProfile::RepeatState;
            for (const State& state : states[period][periods]) {
                const auto max_indel_length = std::max({state.polymorphism_counts.size(), state.error_counts.size(), std::size_t {2}}) - 1;
                for (std::size_t indel_length {1}; indel_length <= max_indel_length; ++indel_length) {
                    os << '\n' << period << ',' << periods << ',' << state.motif << ','
                       << state.reference_count << ',' << state.span << ',' << indel_length << ',';
                    if (indel_length < state.polymorphism_counts.size()) {
                        os << state.polymorphism_counts[indel_length];
                    } else {
                        os << 0;
                    }
                    os << ',';
                    if (indel_length < state.error_counts.size()) {
                        os << state.error_counts[indel_length];
                    } else {
                        os << 0;
                    }
                    os << ',' << state.read_count;
                }
            }
        }
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const IndelProfiler::IndelProfile& indel_profile)
{
    os << indel_profile.states;
    return os;
}

} // namespace octopus
