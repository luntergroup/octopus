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
                                                         boost::optional<ProgressMeter&> progress)
: VariantCallFilter {std::move(facet_factory), measures, std::move(output_config)}
, progress_ {progress}
, annotate_measures_ {output_config.annotate_measures}
{}

void SinglePassVariantCallFilter::filter(const VcfReader& source, VcfWriter& dest, const SampleList& samples) const
{
    assert(dest.is_header_written());
    if (progress_) progress_->start();
    if (can_measure_single_call()) {
        auto p = source.iterate();
        std::for_each(std::move(p.first), std::move(p.second), [&] (const VcfRecord& call) { filter(call, dest); });
    } else {
        for (auto p = source.iterate(); p.first != p.second;) {
            filter(get_next_block(p.first, p.second, samples), dest);
        }
    }
    if (progress_) progress_->stop();
}

void SinglePassVariantCallFilter::filter(const VcfRecord& call, VcfWriter& dest) const
{
    filter(call, measure(call), dest);
}

void SinglePassVariantCallFilter::filter(const std::vector<VcfRecord>& calls, VcfWriter& dest) const
{
    const auto measures = measure(calls);
    assert(measures.size() == calls.size());
    for (auto tup : boost::combine(calls, measures)) {
        filter(tup.get<0>(), tup.get<1>(), dest);
    }
}

template <typename T>
struct MeasureValueVisitor : public boost::static_visitor<T>
{
    template <typename _> T operator()(const boost::optional<_>& value) const { return static_cast<T>(*value); }
    template <typename _> T operator()(const _& value) const noexcept { return static_cast<T>(value); }
};

template <typename T>
auto get_value(const Measure::ResultType& value)
{
    return boost::apply_visitor(MeasureValueVisitor<T> {}, value);
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
