// Copyright (c) 2017 Daniel Cooke
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

std::unique_ptr<Measure> GenotypeQuality::do_clone() const
{
    return std::make_unique<GenotypeQuality>(*this);
}

Measure::ResultType GenotypeQuality::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    std::vector<boost::optional<double>> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        static const std::string gq_field {vcfspec::format::conditionalQuality};
        boost::optional<double> sample_gq {};
        if (call.has_format(gq_field)) {
            sample_gq = std::stod(call.get_sample_value(sample, gq_field).front());
        }
        result.push_back(sample_gq);
    }
    return result;
}

Measure::ResultCardinality GenotypeQuality::do_cardinality() const noexcept
{
    return ResultCardinality::num_samples;
}

std::string GenotypeQuality::do_name() const
{
    return "GQ";
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
