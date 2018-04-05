// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "double_pass_variant_call_filter.hpp"

#include <utility>
#include <iterator>
#include <algorithm>
#include <cassert>

#include "io/variant/vcf_reader.hpp"
#include "io/variant/vcf_writer.hpp"
#include "utils/append.hpp"

namespace octopus { namespace csr {

DoublePassVariantCallFilter::DoublePassVariantCallFilter(FacetFactory facet_factory,
                                                         std::vector<MeasureWrapper> measures,
                                                         OutputOptions output_config,
                                                         ConcurrencyPolicy threading,
                                                         boost::optional<ProgressMeter&> progress)
: VariantCallFilter {std::move(facet_factory), std::move(measures), std::move(output_config), threading}
, info_log_ {logging::InfoLogger {}}
, progress_ {progress}
, current_contig_ {}
{}

void DoublePassVariantCallFilter::filter(const VcfReader& source, VcfWriter& dest, const SampleList& samples) const
{
    assert(dest.is_header_written());
    make_registration_pass(source, samples);
    prepare_for_classification(info_log_);
    make_filter_pass(source, samples, dest);
}

void DoublePassVariantCallFilter::log_registration_pass(Log& log) const
{
    log << "CSR: Starting registration pass";
}

void DoublePassVariantCallFilter::make_registration_pass(const VcfReader& source, const SampleList& samples) const
{
    if (info_log_) log_registration_pass(*info_log_);
    prepare_for_registration(samples);
    if (progress_) progress_->start();
    if (can_measure_multiple_blocks()) {
        std::size_t idx {0};
        for (auto p = source.iterate(); p.first != p.second;) {
            for (auto block : read_next_blocks(p.first, p.second, samples)) {
                record(block, idx, samples);
                idx += block.size();
            }
        }
    } else if (can_measure_single_call()) {
        auto p = source.iterate();
        std::size_t idx {0};
        std::for_each(std::move(p.first), std::move(p.second), [&] (const VcfRecord& call) { record(call, idx++, samples); });
    } else {
        std::size_t idx {0};
        for (auto p = source.iterate(); p.first != p.second;) {
            const auto calls = read_next_block(p.first, p.second, samples);
            record(calls, idx, samples);
            idx += calls.size();
        }
    }
    if (progress_) progress_->stop();
}

void DoublePassVariantCallFilter::record(const VcfRecord& call, const std::size_t record_idx, const SampleList& samples) const
{
    const auto measures_values = measure(call);
    for (std::size_t sample_idx {0}; sample_idx < samples.size(); ++sample_idx) {
        this->record(record_idx, sample_idx, get_sample_values(measures_values, measures_, sample_idx));
    }
    log_progress(mapped_region(call));
}

void DoublePassVariantCallFilter::record(const std::vector<VcfRecord>& calls, std::size_t record_idx, const SampleList& samples) const
{
    if (!calls.empty()) {
        const auto measure_values = measure(calls);
        assert(measure_values.size() == calls.size());
        for (const auto& values : measure_values) {
            for (std::size_t sample_idx {0}; sample_idx < samples.size(); ++sample_idx) {
                this->record(record_idx, sample_idx, get_sample_values(values, measures_, sample_idx));
            }
            ++record_idx;
        }
        log_progress(encompassing_region(calls));
    }
}

void DoublePassVariantCallFilter::log_filter_pass_start(Log& log) const
{
    log << "CSR: Starting filtering pass";
}

namespace {

template <typename Range, typename BinaryPredicate>
bool all_equal(const Range& values, BinaryPredicate pred)
{
    const auto not_pred = [&] (const auto& lhs, const auto& rhs) { return !pred(lhs, rhs); };
    return std::adjacent_find(std::cbegin(values), std::cend(values), not_pred) == std::cend(values);
}

template <typename Range, typename UnaryPredicate>
bool any_of(const Range& values, UnaryPredicate pred)
{
    return std::any_of(std::cbegin(values), std::cend(values), pred);
}

} // namespace

VariantCallFilter::Classification DoublePassVariantCallFilter::merge(const std::vector<Classification>& sample_classifications) const
{
    assert(!sample_classifications.empty());
    if (sample_classifications.size() == 1) {
        return sample_classifications.front();
    }
    Classification result {};
    if (all_equal(sample_classifications, [] (const auto& lhs, const auto& rhs) { return lhs.category == rhs.category; })) {
        result.category = sample_classifications.front().category;
    } else if (any_of(sample_classifications, [] (const auto& c) { return c.category == Classification::Category::unfiltered; })) {
        result.category = Classification::Category::unfiltered;
    } else {
        result.category = Classification::Category::soft_filtered;
    }
    if (result.category != Classification::Category::unfiltered) {
        for (const auto& sample_classification : sample_classifications) {
            utils::append(sample_classification.reasons, result.reasons);
        }
        std::sort(std::begin(result.reasons), std::end(result.reasons));
        result.reasons.erase(std::unique(std::begin(result.reasons), std::end(result.reasons)), std::end(result.reasons));
        result.reasons.shrink_to_fit();
    }
    for (const auto& sample_classification : sample_classifications) {
        if (sample_classification.quality) {
            if (result.quality) {
                result.quality = std::max(*result.quality, *sample_classification.quality);
            } else {
                result.quality = sample_classification.quality;
            }
        }
    }
    return result;
}

void DoublePassVariantCallFilter::make_filter_pass(const VcfReader& source, const SampleList& samples, VcfWriter& dest) const
{
    if (info_log_) log_filter_pass_start(*info_log_);
    if (progress_) {
        progress_->reset();
        progress_->set_max_tick_size(10);
        progress_->start();
    }
    auto p = source.iterate();
    std::size_t idx {0};
    std::for_each(std::move(p.first), std::move(p.second), [&] (const VcfRecord& call) { filter(call, idx++, samples, dest); });
    if (progress_) progress_->stop();
}

std::vector<VariantCallFilter::Classification>
DoublePassVariantCallFilter::classify(std::size_t call_idx, const SampleList& samples) const
{
    std::vector<Classification> result(samples.size());
    for (std::size_t sample_idx {0}; sample_idx < samples.size(); ++sample_idx) {
        result[sample_idx] = this->classify(call_idx, sample_idx);
    }
    return result;
}

void DoublePassVariantCallFilter::filter(const VcfRecord& call, const std::size_t call_idx, const SampleList& samples,
                                         VcfWriter& dest) const
{
    const auto sample_classifications = classify(call_idx, samples);
    const auto call_classification = merge(sample_classifications);
    write(call, call_classification, samples, sample_classifications, dest);
    log_progress(mapped_region(call));
}

static auto expand_lhs_to_zero(const GenomicRegion& region)
{
    return GenomicRegion {region.contig_name(), 0, region.end()};
}

void DoublePassVariantCallFilter::log_progress(const GenomicRegion& region) const
{
    if (progress_) {
        if (current_contig_) {
            if (*current_contig_ != region.contig_name()) {
                progress_->log_completed(*current_contig_);
                current_contig_ = region.contig_name();
            }
        } else {
            current_contig_ = region.contig_name();
        }
        progress_->log_completed(expand_lhs_to_zero(region));
    }
}

} // namespace csr
} // namespace octopus
