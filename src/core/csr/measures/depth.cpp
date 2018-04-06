// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "depth.hpp"

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "../facets/samples.hpp"
#include "../facets/overlapping_reads.hpp"

namespace octopus { namespace csr {

Depth::Depth() : Depth {false, false} {}

Depth::Depth(bool recalculate, bool aggregate_samples)
: recalculate_ {recalculate}
, aggregate_ {aggregate_samples}
{}

std::unique_ptr<Measure> Depth::do_clone() const
{
    return std::make_unique<Depth>(*this);
}

Measure::ResultType Depth::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    if (aggregate_) {
        if (recalculate_) {
            const auto& reads = get_value<OverlappingReads>(facets.at("OverlappingReads"));
            return static_cast<std::size_t>(count_overlapped(reads, call));
        } else {
            return static_cast<std::size_t>(std::stoull(call.info_value(vcfspec::info::combinedReadDepth).front()));
        }
    } else {
        const auto& samples = get_value<Samples>(facets.at("Samples"));
        std::vector<std::size_t> result {};
        result.reserve(samples.size());
        if (recalculate_) {
            const auto& reads = get_value<OverlappingReads>(facets.at("OverlappingReads"));
            for (const auto& sample : samples) {
                result.push_back(count_overlapped(reads.at(sample), call));
            }
        } else {
            for (const auto& sample : samples) {
                result.push_back(std::stoull(call.get_sample_value(sample, vcfspec::format::combinedReadDepth).front()));
            }
        }
        return result;
    }
}

Measure::ResultCardinality Depth::do_cardinality() const noexcept
{
    if (aggregate_) {
        return ResultCardinality::one;
    } else {
        return ResultCardinality::num_samples;
    }
}

std::string Depth::do_name() const
{
    return "DP";
}

std::string Depth::do_describe() const
{
    return "Number of read overlapping the call";
}

std::vector<std::string> Depth::do_requirements() const
{
    std::vector<std::string> result {};
    if (!aggregate_) result.push_back("Samples");
    if (recalculate_) result.push_back("OverlappingReads");
    return result;
}

bool Depth::is_equal(const Measure& other) const noexcept
{
    const auto& other_depth = static_cast<const Depth&>(other);
    return recalculate_ == other_depth.recalculate_ && aggregate_ == other_depth.aggregate_;
}

} // namespace csr
} // namespace octopus
