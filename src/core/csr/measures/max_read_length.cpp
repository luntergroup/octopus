// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "max_read_length.hpp"

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "utils/read_stats.hpp"
#include "../facets/samples.hpp"
#include "../facets/overlapping_reads.hpp"

namespace octopus { namespace csr {

const std::string MaxReadLength::name_ = "MRL";

std::unique_ptr<Measure> MaxReadLength::do_clone() const
{
    return std::make_unique<MaxReadLength>(*this);
}

Measure::ValueType MaxReadLength::get_value_type() const
{
    return std::size_t {};
}

Measure::ResultType MaxReadLength::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& reads = get_value<OverlappingReads>(facets.at("OverlappingReads"));
    Array<ValueType> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        result.emplace_back(max_read_length(reads.at(sample), mapped_region(call)));
    }
    return result;
}

Measure::ResultCardinality MaxReadLength::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& MaxReadLength::do_name() const
{
    return name_;
}

std::string MaxReadLength::do_describe() const
{
    return "Maximum read length overlapping the site";
}

std::vector<std::string> MaxReadLength::do_requirements() const
{
    return {"Samples", "OverlappingReads"};
}

} // namespace csr
} // namespace octopus
