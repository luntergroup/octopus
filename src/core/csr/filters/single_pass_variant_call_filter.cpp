// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "single_pass_variant_call_filter.hpp"

#include <functional>
#include <utility>
#include <iterator>
#include <algorithm>
#include <cassert>

#include <boost/range/combine.hpp>

#include "io/variant/vcf_reader.hpp"
#include "io/variant/vcf_writer.hpp"
#include "utils/append.hpp"

namespace octopus { namespace csr {

SinglePassVariantCallFilter::SinglePassVariantCallFilter(FacetFactory facet_factory,
                                                         std::vector<MeasureWrapper> measures,
                                                         OutputOptions output_config,
                                                         ConcurrencyPolicy threading,
                                                         boost::optional<ProgressMeter&> progress)
: VariantCallFilter {std::move(facet_factory), measures, std::move(output_config), threading}
, progress_ {progress}
, annotate_measures_ {output_config.annotate_measures}
{}

namespace {

template <typename Range, typename BinaryPredicate>
bool all_equal(const Range& values, BinaryPredicate pred)
{
    const auto not_pred = [&](const auto& lhs, const auto& rhs) { return !pred(lhs, rhs); };
    return std::adjacent_find(std::cbegin(values), std::cend(values), not_pred) == std::cend(values);
}

template <typename Range, typename UnaryPredicate>
bool any_of(const Range& values, UnaryPredicate pred)
{
    return std::any_of(std::cbegin(values), std::cend(values), pred);
}

} // namespace

VariantCallFilter::Classification SinglePassVariantCallFilter::classify(const std::vector<Classification>& sample_classifications) const
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

void SinglePassVariantCallFilter::filter(const VcfReader& source, VcfWriter& dest, const SampleList& samples) const
{
    assert(dest.is_header_written());
    if (progress_) progress_->start();
    if (can_measure_multiple_blocks()) {
        for (auto p = source.iterate(); p.first != p.second;) {
            filter(read_next_blocks(p.first, p.second, samples), dest, samples);
        }
    } else if (can_measure_single_call()) {
        auto p = source.iterate();
        std::for_each(std::move(p.first), std::move(p.second), [&] (const VcfRecord& call) { filter(call, dest, samples); });
    } else {
        for (auto p = source.iterate(); p.first != p.second;) {
            filter(read_next_block(p.first, p.second, samples), dest, samples);
        }
    }
    if (progress_) progress_->stop();
}

void SinglePassVariantCallFilter::filter(const VcfRecord& call, VcfWriter& dest, const SampleList& samples) const
{
    filter(call, measure(call), dest, samples);
}

void SinglePassVariantCallFilter::filter(const CallBlock& block, VcfWriter& dest, const SampleList& samples) const
{
    filter(block, measure(block), dest, samples);
}

void SinglePassVariantCallFilter::filter(const std::vector<CallBlock>& blocks, VcfWriter& dest, const SampleList& samples) const
{
    const auto measures = measure(blocks);
    assert(measures.size() == blocks.size());
    for (auto tup : boost::combine(blocks, measures)) {
        filter(tup.get<0>(), tup.get<1>(), dest, samples);
    }
}

void SinglePassVariantCallFilter::filter(const CallBlock& block, const MeasureBlock& measures, VcfWriter& dest, const SampleList& samples) const
{
    assert(measures.size() == block.size());
    for (auto tup : boost::combine(block, measures)) {
        filter(tup.get<0>(), tup.get<1>(), dest, samples);
    }
}

void SinglePassVariantCallFilter::filter(const VcfRecord& call, const MeasureVector& measures, VcfWriter& dest, const SampleList& samples) const
{
    const auto sample_classifications = classify(measures, samples);
    const auto call_classification = classify(sample_classifications);
    if (annotate_measures_) {
        auto annotation_builder = VcfRecord::Builder {call};
        annotate(annotation_builder, measures);
        const auto annotated_call = annotation_builder.build_once();
        write(annotated_call, call_classification, samples, sample_classifications, dest);
    } else {
        write(call, call_classification, samples, sample_classifications, dest);
    }
    log_progress(mapped_region(call));
}

std::vector<VariantCallFilter::Classification>
SinglePassVariantCallFilter::classify(const MeasureVector& call_measures, const SampleList& samples) const
{
    std::vector<Classification> result(samples.size());
    for (std::size_t sample_idx {0}; sample_idx < samples.size(); ++sample_idx) {
        result[sample_idx] = this->classify(get_sample_values(call_measures, measures_, sample_idx));
    }
    return result;
}

static auto expand_lhs_to_zero(const GenomicRegion& region)
{
    return GenomicRegion {region.contig_name(), 0, region.end()};
}

void SinglePassVariantCallFilter::log_progress(const GenomicRegion& region) const
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
