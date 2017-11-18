// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "variant_call_filter.hpp"

#include <utility>
#include <unordered_map>
#include <map>
#include <limits>
#include <numeric>
#include <cmath>

#include <boost/math/distributions/hypergeometric.hpp>

#include "config/common.hpp"
#include "containers/mappable_flat_set.hpp"
#include "containers/mappable_map.hpp"
#include "core/types/variant.hpp"
#include "core/models/haplotype_likelihood_cache.hpp"
#include "basics/aligned_read.hpp"
#include "utils/read_stats.hpp"
#include "utils/maths.hpp"
#include "utils/string_utils.hpp"
#include "utils/genotype_reader.hpp"
#include "utils/append.hpp"
#include "io/variant/vcf_writer.hpp"

namespace octopus { namespace csr  {

auto get_facets(const std::vector<MeasureWrapper>& measures)
{
    std::vector<std::string> result {};
    result.reserve(measures.size()); // Just a guess
    for (const auto& measure : measures) {
        utils::append(measure.requirements(), result);
    }
    std::sort(std::begin(result), std::end(result));
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

// public methods

VariantCallFilter::VariantCallFilter(FacetFactory facet_factory,
                                     std::vector<MeasureWrapper> measures,
                                     OutputOptions output_config)
: facet_factory_ {std::move(facet_factory)}
, facets_ {get_facets(measures)}
, measures_ {std::move(measures)}
, output_config_ {output_config}
{}

void VariantCallFilter::filter(const VcfReader& source, VcfWriter& dest) const
{
    if (!dest.is_header_written()) {
        dest << make_header(source);
    }
    const auto samples = source.fetch_header().samples();
    filter(source, dest, samples);
}

// protected methods

bool VariantCallFilter::can_measure_single_call() const noexcept
{
    return facets_.empty();
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

std::vector<VcfRecord> VariantCallFilter::get_next_block(VcfIterator& first, const VcfIterator& last, const SampleList& samples) const
{
    std::vector<std::pair<VcfRecord, GenomicRegion>> block {};
    for (; first != last; ++first) {
        const VcfRecord& call {*first};
        auto call_phase_region = get_phase_region(call, samples);
        if (!block.empty() && !overlaps(block.back().second, call_phase_region)) {
            return copy_each_first(block);
        }
        block.push_back({call, std::move(call_phase_region)});
    }
    return copy_each_first(block);
}

VariantCallFilter::MeasureVector VariantCallFilter::measure(const VcfRecord& call) const
{
    MeasureVector result(measures_.size());
    std::transform(std::cbegin(measures_), std::cend(measures_), std::begin(result),
                   [&call] (const MeasureWrapper& f) { return f(call); });
    return result;
}

std::vector<VariantCallFilter::MeasureVector> VariantCallFilter::measure(const std::vector<VcfRecord>& calls) const
{
    const auto facets = compute_facets(calls);
    std::vector<MeasureVector> result {};
    result.reserve(calls.size());
    std::transform(std::cbegin(calls), std::cend(calls), std::back_inserter(result),
                   [this, &facets] (const auto& call) { return measure(call, facets); });
    return result;
}

void VariantCallFilter::write(const VcfRecord& call, const Classification& classification, VcfWriter& dest) const
{
    if (classification.category != Classification::Category::hard_filtered) {
        auto filtered_call = construct_template(call);
        annotate(filtered_call, classification);
        dest << filtered_call.build_once();
    }
}

// private methods

VcfHeader VariantCallFilter::make_header(const VcfReader& source) const
{
    VcfHeader::Builder builder {source.fetch_header()};
    if (output_config_.emit_sites_only) {
        builder.clear_format();
    }
    annotate(builder);
    return builder.build_once();
}

VcfRecord::Builder VariantCallFilter::construct_template(const VcfRecord& call) const
{
    VcfRecord::Builder result {call};
    if (output_config_.emit_sites_only) {
        result.clear_format();
    }
    return result;
}

void VariantCallFilter::annotate(VcfRecord::Builder& call, const Classification status) const
{
    if (status.category == Classification::Category::unfiltered) {
        pass(call);
    } else {
        fail(call, std::move(status.reasons));
    }
}

Measure::FacetMap VariantCallFilter::compute_facets(const std::vector<VcfRecord>& calls) const
{
    Measure::FacetMap result {};
    result.reserve(facets_.size());
    for (const auto& facet : facets_) {
        result.emplace(facet, facet_factory_.make(facet, calls));
    }
    return result;
}

VariantCallFilter::MeasureVector VariantCallFilter::measure(const VcfRecord& call, const Measure::FacetMap& facets) const
{
    MeasureVector result(measures_.size());
    std::transform(std::cbegin(measures_), std::cend(measures_), std::begin(result),
                   [&] (const MeasureWrapper& measure) { return measure(call, facets); });
    return result;
}

void VariantCallFilter::pass(VcfRecord::Builder& call) const
{
    call.set_passed();
}

void VariantCallFilter::fail(VcfRecord::Builder& call, std::vector<std::string> reasons) const
{
    for (auto& reason : reasons) {
        call.add_filter(std::move(reason));
    }
}

} // namespace csr
} // namespace octopus
