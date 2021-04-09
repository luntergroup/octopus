// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "variant_call_filter.hpp"

#include <utility>
#include <unordered_map>
#include <map>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <limits>
#include <cmath>
#include <thread>

#include <boost/range/combine.hpp>
#include <boost/multiprecision/gmp.hpp>

#include "config/common.hpp"
#include "containers/mappable_flat_set.hpp"
#include "containers/mappable_map.hpp"
#include "core/types/variant.hpp"
#include "core/models/haplotype_likelihood_array.hpp"
#include "basics/aligned_read.hpp"
#include "utils/read_stats.hpp"
#include "utils/maths.hpp"
#include "utils/string_utils.hpp"
#include "utils/genotype_reader.hpp"
#include "utils/append.hpp"
#include "utils/parallel_transform.hpp"
#include "io/variant/vcf_writer.hpp"
#include "io/variant/vcf_spec.hpp"
#include "core/csr/measures/measure_factory.hpp"

namespace octopus { namespace csr {

namespace {

unsigned get_pool_size(VariantCallFilter::ConcurrencyPolicy policy)
{
    const auto num_cores = std::thread::hardware_concurrency();
    if (policy.max_threads) {
        if (*policy.max_threads > 1) {
            return num_cores > 0 ? std::min(*policy.max_threads, num_cores) : *policy.max_threads;
        } else {
            return 0;
        }
    } else {
        return num_cores > 0 ? num_cores : 8;
    }
}

bool contains(const std::vector<MeasureWrapper>& measures, const std::string& name)
{
    return std::find_if(std::cbegin(measures), std::cend(measures),
                        [&] (const auto& measure) { return measure.name() == name; }) != std::cend(measures);
}

} // namespace

// public methods

VariantCallFilter::VariantCallFilter(FacetFactory facet_factory,
                                     std::vector<MeasureWrapper> measures,
                                     OutputOptions output_config,
                                     ConcurrencyPolicy threading)
: measures_ {std::move(measures)}
, debug_log_ {logging::get_debug_log()}
, facet_factory_ {std::move(facet_factory)}
, facet_names_ {}
, output_config_ {output_config}
, duplicate_measures_ {}
, workers_ {get_pool_size(threading)}
{
    std::unordered_map<MeasureWrapper, int> measure_counts {};
    measure_counts.reserve(measures_.size());
    for (const auto& m : measures_) {
        ++measure_counts[m];
        if (measure_counts[m] == 2) {
            duplicate_measures_.push_back(m);
        }
    }
    duplicate_measures_.shrink_to_fit();
    logging::WarningLogger warn_log {};
    for (const auto& annotation : output_config_.annotations) {
        if (!contains(measures_, annotation)) {
            measures_.push_back(make_measure(annotation));
        }
    }
    measures_.shrink_to_fit();
    facet_names_ = get_all_requirements(measures_);
}

inline VcfHeader read_header(const VcfWriter& writer)
{
    assert(writer.is_header_written());
    assert(writer.path());
    VcfReader tmp {*writer.path()};
    return tmp.fetch_header();
}

std::string VariantCallFilter::name() const
{
    return do_name();
}

void VariantCallFilter::filter(const VcfReader& source, VcfWriter& dest, boost::optional<VcfHeader> template_header) const
{
    if (dest.is_header_written()) {
        dest.close();
        const auto header = make_header(read_header(dest));
        dest.open(true);
        dest << header;
        filter(source, dest, header);
    } else if (template_header) {
        const auto header = make_header(*template_header);
        dest << header;
        filter(source, dest, header);
    } else {
        const auto header = make_header(source);
        dest << header;
        filter(source, dest, header);
    }
}

// protected methods

namespace {
template <typename Range, typename BinaryPredicate>
bool all_equal(const Range& values, BinaryPredicate pred)
{
    const auto not_pred = [&](const auto& lhs, const auto& rhs) { return !pred(lhs, rhs); };
    return std::adjacent_find(std::cbegin(values), std::cend(values), not_pred) == std::cend(values);
}
} // namespace

VariantCallFilter::Classification
VariantCallFilter::merge(const ClassificationList& sample_classifications, const MeasureVector& measures) const
{
    assert(!sample_classifications.empty());
    if (sample_classifications.size() == 1) {
        return sample_classifications.front();
    }
    Classification result {};
    result.quality = compute_joint_quality(sample_classifications, measures);
    if (!result.quality && all_equal(sample_classifications, [] (const auto& lhs, const auto& rhs) { return lhs.category == rhs.category; })) {
        result.category = sample_classifications.front().category;
    } else if (is_soft_filtered(sample_classifications, result.quality, measures, result.reasons)) {
        result.category = Classification::Category::soft_filtered;
    } else {
        result.category = Classification::Category::unfiltered;
    }
    if (result.category != Classification::Category::unfiltered && result.reasons.empty()) {
        result.reasons = compute_reason_union(sample_classifications);
    }
    return result;
}

VariantCallFilter::Classification VariantCallFilter::merge(const ClassificationList& sample_classifications) const
{
    return this->merge(sample_classifications, {});
}

bool VariantCallFilter::can_measure_single_call() const noexcept
{
    return facet_names_.empty();
}

bool VariantCallFilter::can_measure_multiple_blocks() const noexcept
{
    return is_multithreaded();
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

template <typename Range>
bool can_add_to_phase_block(const VcfRecord& call, const GenomicRegion& call_phase_region, const Range& block)
{
    if (block.empty()) return true;
    if (overlaps(block.back().second, call_phase_region)) return true;
    return is_refcall(call) && is_same_contig(call, block.back().first);
}

} // namespace

VariantCallFilter::CallBlock
VariantCallFilter::read_next_block(VcfIterator& first, const VcfIterator& last, const SampleList& samples) const
{
    std::vector<std::pair<VcfRecord, GenomicRegion>> block {};
    for (; first != last; ++first) {
        const VcfRecord& call {*first};
        auto call_phase_region = get_phase_region(call, samples);
        if (!can_add_to_phase_block(call, call_phase_region, block)) {
            return copy_each_first(block);
        }
        block.emplace_back(call, std::move(call_phase_region));
    }
    return copy_each_first(block);
}

std::vector<VariantCallFilter::CallBlock>
VariantCallFilter::read_next_blocks(VcfIterator& first, const VcfIterator& last, const SampleList& samples) const
{
    std::vector<VariantCallFilter::CallBlock> result {};
    if (can_measure_multiple_blocks()) {
        const auto max_blocks = max_concurrent_blocks();
        result.reserve(max_blocks);
        while (result.size() < max_blocks) {
            result.push_back(read_next_block(first, last, samples));
            if (first == last || !is_same_contig(*first, result.back().front())) break;
        }
    } else {
        result.push_back(read_next_block(first, last, samples));
    }
    return result;
}

VariantCallFilter::MeasureVector VariantCallFilter::measure(const VcfRecord& call) const
{
    MeasureVector result(measures_.size());
    if (duplicate_measures_.empty()) {
        std::transform(std::cbegin(measures_), std::cend(measures_), std::begin(result),
                       [&call] (const MeasureWrapper& m) { return m(call); });
    } else {
        std::unordered_map<MeasureWrapper, Measure::ResultType> result_buffer {};
        result_buffer.reserve(duplicate_measures_.size());
        for (const auto& m : duplicate_measures_) {
            result_buffer.emplace(m, m(call));
        }
        std::transform(std::cbegin(measures_), std::cend(measures_), std::begin(result),
                       [&call, &result_buffer] (const MeasureWrapper& m) -> Measure::ResultType {
                           auto itr = result_buffer.find(m);
                           if (itr != std::cend(result_buffer)) {
                               return itr->second;
                           } else {
                               return m(call);
                           }
                       });
    }
    return result;
}

VariantCallFilter::MeasureBlock VariantCallFilter::measure(const CallBlock& block) const
{
    const auto facets = compute_facets(block);
    auto result = measure(block, facets);
    return result;
}

std::vector<VariantCallFilter::MeasureBlock> VariantCallFilter::measure(const std::vector<CallBlock>& blocks) const
{
    std::vector<MeasureBlock> result {};
    result.reserve(blocks.size());
    if (is_multithreaded()) {
        const auto facets = compute_facets(blocks);
        if (debug_log_) {
            stream(*debug_log_) << "Measuring " << blocks.size() << " blocks with " << workers_.size() << " threads";
        }
        transform(std::cbegin(blocks), std::cend(blocks), std::cbegin(facets), std::back_inserter(result),
                  [this] (const auto& block, const auto& block_facets) {
                      return this->measure(block, block_facets);
                  }, workers_);
    } else {
        for (const CallBlock& block : blocks) {
            result.push_back(measure(block));
        }
    }
    return result;
}

void VariantCallFilter::write(const VcfRecord& call, const Classification& classification, VcfWriter& dest) const
{
    if (!is_hard_filtered(classification)) {
        auto filtered_call = construct_template(call);
        annotate(filtered_call, classification);
        dest << filtered_call.build_once();
    }
}

void VariantCallFilter::write(const VcfRecord& call, const Classification& classification,
                              const SampleList& samples, const ClassificationList& sample_classifications,
                              VcfWriter& dest) const
{
    if (!is_hard_filtered(classification)) {
        auto filtered_call = construct_template(call);
        annotate(filtered_call, classification);
        annotate(filtered_call, samples, sample_classifications);
        dest << filtered_call.build_once();
    }
}

bool VariantCallFilter::measure_annotations_requested() const noexcept
{
    return output_config_.annotate_all_active_measures || !output_config_.annotations.empty();
}

void VariantCallFilter::annotate(VcfRecord::Builder& call, const MeasureVector& measures, const VcfHeader& header) const
{
    if (output_config_.clear_info) {
        call.clear_info();
    }
    for (auto p : boost::combine(measures_, measures)) {
        const MeasureWrapper& measure {p.get<0>()};
        if (is_requested_annotation(measure)) {
            const Measure::ResultType& measured_value {p.get<1>()};
            measure.annotate(call, measured_value, header, output_config_.aggregate_allele_annotations);
        }
    }
}

// private methods

boost::optional<Phred<double>>
VariantCallFilter::compute_joint_quality(const ClassificationList& sample_classifications, const MeasureVector& measures) const
{
    std::vector<Phred<double>> sample_classification_qualities {};
    sample_classification_qualities.reserve(sample_classifications.size());
    for (const auto& sample_classification : sample_classifications) {
        if (sample_classification.quality) {
            sample_classification_qualities.push_back(*sample_classification.quality);
        }
    }
    if (!sample_classification_qualities.empty()) {
        return compute_joint_quality(sample_classification_qualities);
    }
    return boost::none;
}

bool
VariantCallFilter::is_soft_filtered(const ClassificationList& sample_classifications,
                                    boost::optional<Phred<double>> joint_quality,
                                    const MeasureVector& measures,
                                    std::vector<std::string>& reasons) const
{
    return std::all_of(std::cbegin(sample_classifications), std::cend(sample_classifications),
                       [] (const auto& c) { return c.category != Classification::Category::unfiltered; });
}

Phred<double> ln_probability_true_to_phred(const double ln_prob_true)
{
    assert(ln_prob_true < 0);
    if (ln_prob_true <= -1e-10) {
        return Phred<double> {maths::probability_true_to_phred<double>(std::exp(ln_prob_true))};
    } else {
        using BigFloat = boost::multiprecision::mpf_float_500;
        return octopus::ln_probability_true_to_phred<double>(BigFloat {ln_prob_true});
    }
}

Phred<double> VariantCallFilter::compute_joint_quality(const std::vector<Phred<double>>& qualities) const
{
    assert(!qualities.empty());
    if (qualities.size() == 1) return qualities.front();
    if (std::any_of(std::cbegin(qualities), std::cend(qualities), [] (auto p) { return p.score() <= 0; })) {
        return Phred<double> {0.0};
    }
    std::vector<double> log_probs(qualities.size());
    std::transform(std::cbegin(qualities), std::cend(qualities), std::begin(log_probs),
                   [] (auto p) { return std::log(p.probability_true()); });
    const auto ln_prob_all_good = std::accumulate(std::cbegin(log_probs), std::cend(log_probs), 0.0);
    return ln_probability_true_to_phred(ln_prob_all_good / log_probs.size());
}

std::vector<std::string> VariantCallFilter::compute_reason_union(const ClassificationList& sample_classifications) const
{
    std::vector<std::string> result {};
    for (const auto& sample_classification : sample_classifications) {
        utils::append(sample_classification.reasons, result);
    }
    std::sort(std::begin(result), std::end(result));
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

VcfHeader VariantCallFilter::make_header(const VcfHeader& source) const
{
    VcfHeader::Builder builder {source};
    if (output_config_.clear_info) {
        builder.clear_info();
    }
    if (measure_annotations_requested()) {
        for (const auto& measure : measures_) {
            if (is_requested_annotation(measure)) {
                measure.annotate(builder, output_config_.aggregate_allele_annotations);
            }
        }
    }
    if (output_config_.emit_sites_only) {
        builder.clear_format();
    }
    annotate(builder);
    return builder.build_once();
}

VcfHeader VariantCallFilter::make_header(const VcfReader& source) const
{
    return make_header(source.fetch_header());
}

VcfRecord::Builder VariantCallFilter::construct_template(const VcfRecord& call) const
{
    VcfRecord::Builder result {call};
    if (output_config_.emit_sites_only) {
        result.clear_format();
    }
    if (output_config_.clear_existing_filters) {
        result.clear_filter();
        result.clear_all_sample_filters();
    }
    return result;
}

bool VariantCallFilter::is_requested_annotation(const MeasureWrapper& measure) const noexcept
{
    return output_config_.annotate_all_active_measures || output_config_.annotations.count(measure.name()) == 1;
}

bool VariantCallFilter::is_hard_filtered(const Classification& classification) const noexcept
{
    return classification.category == Classification::Category::hard_filtered;
}

void VariantCallFilter::annotate(VcfRecord::Builder& call, const SampleList& samples, const ClassificationList& sample_classifications) const
{
    assert(samples.size() == sample_classifications.size());
    bool all_hard_filtered {true};
    auto quality_name = this->genotype_quality_name();
    if (quality_name) {
        call.add_format(std::move(*quality_name));
    }
    for (auto p : boost::combine(samples, sample_classifications)) {
        const SampleName& sample {p.get<0>()};
        const Classification& sample_classification {p.get<1>()};
        if (!is_hard_filtered(sample_classification)) {
            annotate(call, sample, sample_classification);
            all_hard_filtered = false;
        } else {
            call.clear_format(sample);
        }
    }
    if (all_hard_filtered) {
        call.clear_format();
    } else {
        call.add_format(vcfspec::format::filter);
    }
}

void VariantCallFilter::annotate(VcfRecord::Builder& call, const SampleName& sample, Classification status) const
{
    if (status.category == Classification::Category::unfiltered) {
        pass(sample, call);
    } else {
        fail(sample, call, std::move(status.reasons));
    }
    const auto quality_name = this->genotype_quality_name();
    if (quality_name) {
        if (status.quality) {
            call.set_format(sample, *quality_name, utils::to_string(status.quality->score(), 2));
        } else {
            call.set_format_missing(sample, *quality_name);
        }
    }
}

void VariantCallFilter::annotate(VcfRecord::Builder& call, const Classification status) const
{
    if (status.category == Classification::Category::unfiltered) {
        pass(call);
    } else {
        fail(call, std::move(status.reasons));
    }
    auto quality_name = this->call_quality_name();
    if (quality_name) {
        call.add_info(*quality_name);
        if (status.quality) {
            call.set_info(*quality_name, utils::to_string(status.quality->score(), 3, utils::PrecisionRule::sf));
        } else {
            call.set_info_missing(*quality_name);
        }
    }
}

auto make_map(const std::vector<std::string>& names, std::vector<FacetWrapper>&& facets)
{
    assert(names.size() == facets.size());
    Measure::FacetMap result {};
    result.reserve(names.size());
    for (auto tup : boost::combine(names, std::move(facets))) {
        result.emplace(tup.get<0>(), std::move(tup.get<1>()));
    }
    return result;
}

Measure::FacetMap VariantCallFilter::compute_facets(const CallBlock& block) const
{
    return make_map(facet_names_, facet_factory_.make(facet_names_, block));
}

std::vector<Measure::FacetMap> VariantCallFilter::compute_facets(const std::vector<CallBlock>& blocks) const
{
    if (debug_log_ && !blocks.empty()) {
        const auto blocks_region = closed_region(blocks.front().front(), blocks.back().back());
        stream(*debug_log_) << "Computing facets in blocks region " << blocks_region << " containing " << blocks.size() << " blocks";
    }
    auto facets = facet_factory_.make(facet_names_, blocks, workers_);
    std::vector<Measure::FacetMap> result {};
    result.reserve(blocks.size());
    for (auto& block : facets) {
        result.push_back(make_map(facet_names_, std::move(block)));
    }
    return result;
}

VariantCallFilter::MeasureBlock VariantCallFilter::measure(const CallBlock& block, const Measure::FacetMap& facets) const
{
    if (debug_log_ && !block.empty()) {
        stream(*debug_log_) << "Measuring block " << encompassing_region(block) << " containing " << block.size() << " calls";
    }
    MeasureBlock result(block.size());
    std::transform(std::cbegin(block), std::cend(block), std::begin(result),
                   [&] (const VcfRecord& call) { return this->measure(call, facets); });
    return result;
}

VariantCallFilter::MeasureVector VariantCallFilter::measure(const VcfRecord& call, const Measure::FacetMap& facets) const
{
    MeasureVector result(measures_.size());
    if (duplicate_measures_.empty()) {
        std::transform(std::cbegin(measures_), std::cend(measures_), std::begin(result),
                       [&] (const MeasureWrapper& m) { return m(call, facets); });
    } else {
        std::unordered_map<MeasureWrapper, Measure::ResultType> result_buffer {};
        result_buffer.reserve(duplicate_measures_.size());
        for (const auto& m : duplicate_measures_) {
            result_buffer.emplace(m, m(call, facets));
        }
        std::transform(std::cbegin(measures_), std::cend(measures_), std::begin(result),
                       [&] (const MeasureWrapper& m) -> Measure::ResultType {
                           auto itr = result_buffer.find(m);
                           if (itr != std::cend(result_buffer)) {
                               return itr->second;
                           } else {
                               return m(call, facets);
                           }
                       });
    }
    return result;
}

void VariantCallFilter::pass(const SampleName& sample, VcfRecord::Builder& call) const
{
    call.set_passed(sample);
}

void VariantCallFilter::pass(VcfRecord::Builder& call) const
{
    call.set_passed();
}

void VariantCallFilter::fail(const SampleName& sample, VcfRecord::Builder& call, std::vector<std::string> reasons) const
{
    for (auto& reason : reasons) {
        call.add_filter(sample, std::move(reason));
    }
}

void VariantCallFilter::fail(VcfRecord::Builder& call, std::vector<std::string> reasons) const
{
    for (auto& reason : reasons) {
        call.add_filter(std::move(reason));
    }
}

bool VariantCallFilter::is_multithreaded() const noexcept
{
    return !workers_.empty();
}

unsigned VariantCallFilter::max_concurrent_blocks() const noexcept
{
    if (is_multithreaded()) {
        return std::min(100 * workers_.size(), std::size_t {10'000});
    } else {
        return 1;
    }
}

} // namespace csr
} // namespace octopus
