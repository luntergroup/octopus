// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "mapping_quality_zero_count.hpp"

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "utils/read_stats.hpp"
#include "../facets/overlapping_reads.hpp"

namespace octopus { namespace csr {

const std::string MappingQualityZeroCount::name_ = "MQ0";

MappingQualityZeroCount::MappingQualityZeroCount(bool recalculate) : recalculate_ {recalculate} {}

std::unique_ptr<Measure> MappingQualityZeroCount::do_clone() const
{
    return std::make_unique<MappingQualityZeroCount>(*this);
}

Measure::ResultType MappingQualityZeroCount::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    if (recalculate_) {
        const auto& reads = get_value<OverlappingReads>(facets.at("OverlappingReads"));
        return count_mapq_zero(reads);
    } else {
        return static_cast<std::size_t>(std::stoull(call.info_value("MQ0").front()));
    }
}

Measure::ResultCardinality MappingQualityZeroCount::do_cardinality() const noexcept
{
    return ResultCardinality::one;
}

const std::string& MappingQualityZeroCount::do_name() const
{
    return name_;
}

std::string MappingQualityZeroCount::do_describe() const
{
    return "Number of reads overlapping the call with mapping quality zero";
}

std::vector<std::string> MappingQualityZeroCount::do_requirements() const
{
    if (recalculate_) {
        return {"OverlappingReads"};
    } else {
        return {};
    }
}

bool MappingQualityZeroCount::is_equal(const Measure& other) const noexcept
{
    return recalculate_ == static_cast<const MappingQualityZeroCount&>(other).recalculate_;
}

} // namespace csr
} // namespace octopus
