// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "double_pass_variant_call_filter.hpp"

#include <utility>
#include <iterator>
#include <algorithm>
#include <cassert>

#include "io/variant/vcf_reader.hpp"
#include "io/variant/vcf_writer.hpp"

namespace octopus { namespace csr  {

DoublePassVariantCallFilter::DoublePassVariantCallFilter(FacetFactory facet_factory,
                                                         std::vector<MeasureWrapper> measures,
                                                         OutputOptions output_config,
                                                         boost::optional<ProgressMeter&> progress)
: VariantCallFilter {std::move(facet_factory), std::move(measures), std::move(output_config)}
, progress_ {progress}
{}

void DoublePassVariantCallFilter::filter(const VcfReader& source, VcfWriter& dest, const SampleList& samples) const
{
    assert(dest.is_header_written());
    make_registration_pass(source, samples);
    make_filter_pass(source, dest);
}

void DoublePassVariantCallFilter::make_registration_pass(const VcfReader& source, const SampleList& samples) const
{
    if (can_measure_single_call()) {
        auto p = source.iterate();
        std::size_t idx {0};
        std::for_each(std::move(p.first), std::move(p.second), [&] (const VcfRecord& call) { record(call, idx++); });
    } else {
        std::size_t idx {0};
        for (auto p = source.iterate(); p.first != p.second;) {
            const auto calls = get_next_block(p.first, p.second, samples);
            record(calls, idx);
            idx += calls.size();
        }
    }
}

void DoublePassVariantCallFilter::record(const VcfRecord& call, const std::size_t idx) const
{
    record(idx, measure(call));
}

void DoublePassVariantCallFilter::record(const std::vector<VcfRecord>& calls, std::size_t first_idx) const
{
    const auto measures = measure(calls);
    assert(measures.size() == calls.size());
    for (const auto& m : measures) {
        record(first_idx++, m);
    }
}

void DoublePassVariantCallFilter::make_filter_pass(const VcfReader& source, VcfWriter& dest) const
{
    auto p = source.iterate();
    std::size_t idx {0};
    std::for_each(std::move(p.first), std::move(p.second), [&] (const VcfRecord& call) { filter(call, idx++, dest); });
}

void DoublePassVariantCallFilter::filter(const VcfRecord& call, const std::size_t idx, VcfWriter& dest) const
{
    write(call, classify(idx), dest);
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
