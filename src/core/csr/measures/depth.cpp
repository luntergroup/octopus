// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "depth.hpp"

#include <boost/variant.hpp>
#include <boost/lexical_cast.hpp>

#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "../facets/samples.hpp"
#include "../facets/overlapping_reads.hpp"

namespace octopus { namespace csr {

const std::string Depth::name_ = "DP";

Depth::Depth() : Depth {false, false} {}

Depth::Depth(bool recalculate, bool aggregate_samples)
: recalculate_ {recalculate}
, aggregate_ {aggregate_samples}
{}

std::unique_ptr<Measure> Depth::do_clone() const
{
    return std::make_unique<Depth>(*this);
}

Measure::ValueType Depth::get_value_type() const
{
    return std::size_t {};
}

Measure::ResultType Depth::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    if (aggregate_) {
        ValueType result {};
        if (recalculate_) {
            const auto& reads = get_value<OverlappingReads>(facets.at("OverlappingReads"));
            result = static_cast<std::size_t>(count_overlapped(reads, call));
        } else {
            result = boost::lexical_cast<std::size_t>(call.info_value(vcfspec::info::combinedReadDepth).front());
        }
        return result;
    } else {
        const auto& samples = get_value<Samples>(facets.at("Samples"));
        Array<ValueType> result {};
        result.reserve(samples.size());
        if (recalculate_) {
            const auto& reads = get_value<OverlappingReads>(facets.at("OverlappingReads"));
            for (const auto& sample : samples) {
                result.emplace_back(static_cast<std::size_t>(count_overlapped(reads.at(sample), call)));
            }
        } else {
            for (const auto& sample : samples) {
                result.emplace_back(boost::lexical_cast<std::size_t>(call.get_sample_value(sample, vcfspec::format::combinedReadDepth).front()));
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
        return ResultCardinality::samples;
    }
}

const std::string& Depth::do_name() const
{
    return name_;
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
