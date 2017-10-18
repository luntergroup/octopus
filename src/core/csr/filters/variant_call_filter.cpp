// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "variant_call_filter.hpp"

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
#include "basics/genomic_region.hpp"
#include "basics/aligned_read.hpp"
#include "utils/read_stats.hpp"
#include "utils/maths.hpp"
#include "utils/string_utils.hpp"
#include "utils/append.hpp"
#include "io/variant/vcf_reader.hpp"
#include "io/variant/vcf_writer.hpp"
#include "io/variant/vcf_header.hpp"
#include "../utils/genotype_reader.hpp"

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

VariantCallFilter::VariantCallFilter(const ReferenceGenome& reference, std::vector<MeasureWrapper> measures)
: reference_ {reference}
, measures_ {std::move(measures)}
, facets_ {get_facets(measures_)}
{}

namespace {

GenomicRegion get_phase_set(const VcfRecord& record, const SampleName& sample)
{
    const std::string& begin{record.get_sample_value(sample, "PS").front()};
    return GenomicRegion {record.chrom(), static_cast<GenomicRegion::Position>(std::stoi(begin)), record.pos()};
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
std::vector<T> copy_first(const std::vector<std::pair<T, _>>& items)
{
    std::vector<T> result {};
    result.reserve(items.size());
    std::transform(std::cbegin(items), std::cend(items), std::back_inserter(result),
                   [] (const auto& p) { return p.first; });
    return result;
}

} // namespace

void VariantCallFilter::filter(const VcfReader& source, VcfWriter& dest)
{
    if (!dest.is_header_written()) {
        auto header = source.fetch_header();
        annotate(header);
        dest << header;
    }
    const auto samples = source.fetch_header().samples();
    if (facets_.empty()) {
        auto p = source.iterate();
        std::for_each(std::move(p.first), std::move(p.second), [&] (const VcfRecord& call) {
            dest << filter(call);
        });
    } else {
        std::vector<std::pair<VcfRecord, GenomicRegion>> block {};
        for (auto p = source.iterate(); p.first != p.second; ++p.first) {
            const VcfRecord& call {*p.first};
            auto call_phase_region = get_phase_region(call, samples);
            if (!block.empty() && !overlaps(block.back().second, call_phase_region)) {
                dest << filter(copy_first(block));
                block.clear();
            }
            block.push_back({call, std::move(call_phase_region)});
        }
        dest << filter(copy_first(block));
    }
}

VcfRecord VariantCallFilter::filter(const VcfRecord& call) const
{
    VcfRecord::Builder filtered_call {call};
    annotate(filtered_call, classify(measure(call)));
    return filtered_call.build_once();
}

std::vector<VcfRecord> VariantCallFilter::filter(std::vector<VcfRecord> calls) const
{
    return {};
}

void VariantCallFilter::annotate(VcfRecord::Builder& call, const Classification status) const
{
    if (status.category != Classification::Category::filtered) {
        pass(call);
    } else {
        fail(call, std::move(status.reasons));
    }
}

VariantCallFilter::FacetSet VariantCallFilter::compute_facets(const VcfRecord& call) const
{
    return {};
}

VariantCallFilter::MeasureVector VariantCallFilter::measure(const VcfRecord& call) const
{
    MeasureVector result(measures_.size());
    std::transform(std::cbegin(measures_), std::cend(measures_), std::begin(result),
                   [&call] (const MeasureWrapper& measure) -> MeasureDomain {
                       return measure(call);
                   });
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
