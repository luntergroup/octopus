// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "single_pass_variant_call_filter.hpp"

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

void SinglePassVariantCallFilter::filter(const VcfReader& source, VcfWriter& dest, const SampleList& samples) const
{
    assert(dest.is_header_written());
    if (progress_) progress_->start();
    if (can_measure_multiple_blocks()) {
        for (auto p = source.iterate(); p.first != p.second;) {
            filter(read_next_blocks(p.first, p.second, samples), dest);
        }
    } else if (can_measure_single_call()) {
        auto p = source.iterate();
        std::for_each(std::move(p.first), std::move(p.second), [&] (const VcfRecord& call) { filter(call, dest); });
    } else {
        for (auto p = source.iterate(); p.first != p.second;) {
            filter(read_next_block(p.first, p.second, samples), dest);
        }
    }
    if (progress_) progress_->stop();
}

void SinglePassVariantCallFilter::filter(const VcfRecord& call, VcfWriter& dest) const
{
    filter(call, measure(call), dest);
}

void SinglePassVariantCallFilter::filter(const CallBlock& block, VcfWriter& dest) const
{
    filter(block, measure(block), dest);
}

void SinglePassVariantCallFilter::filter(const std::vector<CallBlock>& blocks, VcfWriter& dest) const
{
    const auto measures = measure(blocks);
    assert(measures.size() == blocks.size());
    for (auto tup : boost::combine(blocks, measures)) {
        filter(tup.get<0>(), tup.get<1>(), dest);
    }
}

void SinglePassVariantCallFilter::filter(const CallBlock& block, const MeasureBlock& measures, VcfWriter& dest) const
{
    assert(measures.size() == block.size());
    for (auto tup : boost::combine(block, measures)) {
        filter(tup.get<0>(), tup.get<1>(), dest);
    }
}

void SinglePassVariantCallFilter::filter(const VcfRecord& call, const MeasureVector& measures, VcfWriter& dest) const
{
    if (annotate_measures_) {
        auto annotation_builder = VcfRecord::Builder {call};
        annotate(annotation_builder, measures);
        const auto annotated_call = annotation_builder.build_once();
        write(annotated_call, classify(measures), dest);
    } else {
        write(call, classify(measures), dest);
    }
    log_progress(mapped_region(call));
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
