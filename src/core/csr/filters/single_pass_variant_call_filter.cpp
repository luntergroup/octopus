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

void SinglePassVariantCallFilter::filter(const VcfReader& source, VcfWriter& dest, const VcfHeader& header) const
{
    assert(dest.is_header_written());
    if (progress_) progress_->start();
    if (can_measure_multiple_blocks()) {
        for (auto p = source.iterate(); p.first != p.second;) {
            filter(read_next_blocks(p.first, p.second, header.samples()), dest, header);
        }
    } else if (can_measure_single_call()) {
        auto p = source.iterate();
        std::for_each(std::move(p.first), std::move(p.second), [&] (const VcfRecord& call) { filter(call, dest, header); });
    } else {
        for (auto p = source.iterate(); p.first != p.second;) {
            filter(read_next_block(p.first, p.second, header.samples()), dest, header);
        }
    }
    if (progress_) progress_->stop();
}

void SinglePassVariantCallFilter::filter(const VcfRecord& call, VcfWriter& dest, const VcfHeader& header) const
{
    filter(call, measure(call), dest, header);
}

void SinglePassVariantCallFilter::filter(const CallBlock& block, VcfWriter& dest, const VcfHeader& header) const
{
    filter(block, measure(block), dest, header);
}

void SinglePassVariantCallFilter::filter(const std::vector<CallBlock>& blocks, VcfWriter& dest, const VcfHeader& header) const
{
    const auto measures = measure(blocks);
    assert(measures.size() == blocks.size());
    for (auto tup : boost::combine(blocks, measures)) {
        filter(tup.get<0>(), tup.get<1>(), dest, header);
    }
}

void SinglePassVariantCallFilter::filter(const CallBlock& block, const MeasureBlock& measures, VcfWriter& dest, const VcfHeader& header) const
{
    assert(measures.size() == block.size());
    for (auto tup : boost::combine(block, measures)) {
        filter(tup.get<0>(), tup.get<1>(), dest, header);
    }
}

void SinglePassVariantCallFilter::filter(const VcfRecord& call, const MeasureVector& measures, VcfWriter& dest, const VcfHeader& header) const
{
    const auto sample_classifications = classify(measures, header.samples());
    const auto call_classification = merge(sample_classifications, measures);
    if (annotate_measures_) {
        VcfRecord::Builder annotation_builder {call};
        annotate(annotation_builder, measures, header);
        write(annotation_builder.build_once(), call_classification, header.samples(), sample_classifications, dest);
    } else {
        write(call, call_classification, header.samples(), sample_classifications, dest);
    }
    log_progress(mapped_region(call));
}

VariantCallFilter::ClassificationList
SinglePassVariantCallFilter::classify(const MeasureVector& call_measures, const SampleList& samples) const
{
    ClassificationList result(samples.size());
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
