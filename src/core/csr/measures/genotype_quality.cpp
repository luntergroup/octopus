// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "genotype_quality.hpp"

#include <utility>
#include <algorithm>

#include <boost/variant.hpp>
#include <boost/optional.hpp>

#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "../facets/samples.hpp"

namespace octopus { namespace csr {

const std::string GenotypeQuality::name_ = "GQ";

std::unique_ptr<Measure> GenotypeQuality::do_clone() const
{
    return std::make_unique<GenotypeQuality>(*this);
}

Measure::ValueType GenotypeQuality::get_value_type() const
{
    return double {};
}

Measure::ResultType GenotypeQuality::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    Array<Optional<ValueType>> result(samples.size());
    for (std::size_t s {0}; s < samples.size(); ++s) {
        static const std::string gq_field {vcfspec::format::conditionalQuality};
        if (call.has_format(gq_field)) {
            result[s] = std::stod(call.get_sample_value(samples[s], gq_field).front());
        }
    }
    return result;
}

Measure::ResultCardinality GenotypeQuality::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& GenotypeQuality::do_name() const
{
    return name_;
}

std::string GenotypeQuality::do_describe() const
{
    return "GQ of each sample";
}

std::vector<std::string> GenotypeQuality::do_requirements() const
{
    return {"Samples"};
}

} // namespace csr
} // namespace octopus
