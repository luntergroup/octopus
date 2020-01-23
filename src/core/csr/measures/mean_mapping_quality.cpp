// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "mean_mapping_quality.hpp"

#include <cassert>

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "utils/read_stats.hpp"
#include "io/variant/vcf_spec.hpp"
#include "../facets/overlapping_reads.hpp"

namespace octopus { namespace csr {

const std::string MeanMappingQuality::name_ = "MQ";

MeanMappingQuality::MeanMappingQuality(bool recalculate) : recalculate_ {recalculate} {}

std::unique_ptr<Measure> MeanMappingQuality::do_clone() const
{
    return std::make_unique<MeanMappingQuality>(*this);
}

Measure::ResultType MeanMappingQuality::get_default_result() const
{
    return double {};
}

Measure::ResultType MeanMappingQuality::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    if (recalculate_) {
        const auto& reads = get_value<OverlappingReads>(facets.at("OverlappingReads"));
        assert(!reads.empty());
        return rmq_mapping_quality(reads, mapped_region(call));
    } else {
        return std::stod(call.info_value(vcfspec::info::rmsMappingQuality).front());
    }
}

Measure::ResultCardinality MeanMappingQuality::do_cardinality() const noexcept
{
    return ResultCardinality::one;
}

const std::string& MeanMappingQuality::do_name() const
{
    return name_;
}

std::string MeanMappingQuality::do_describe() const
{
    return "Mean mapping quality of reads overlapping the call";
}

std::vector<std::string> MeanMappingQuality::do_requirements() const
{
    if (recalculate_) {
        return {"OverlappingReads"};
    } else {
        return {};
    }
}

bool MeanMappingQuality::is_equal(const Measure& other) const noexcept
{
    return recalculate_ == static_cast<const MeanMappingQuality&>(other).recalculate_;
}

} // namespace csr
} // namespace octopus
