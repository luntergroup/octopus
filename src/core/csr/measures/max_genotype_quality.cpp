// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "max_genotype_quality.hpp"

#include <utility>
#include <algorithm>

#include <boost/variant.hpp>
#include <boost/optional.hpp>

#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "../facets/samples.hpp"

namespace octopus { namespace csr {

std::unique_ptr<Measure> MaxGenotypeQuality::do_clone() const
{
    return std::make_unique<MaxGenotypeQuality>(*this);
}

Measure::ResultType MaxGenotypeQuality::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    boost::optional<double> result {};
    for (const auto& sample : samples) {
        static const std::string gq_field {vcfspec::format::conditionalQuality};
        if (call.has_format(gq_field)) {
            const auto sample_gq = std::stod(call.get_sample_value(sample, gq_field).front());
            if (result) {
                result = std::max(sample_gq, *result);
            } else {
                result = sample_gq;
            }
        }
    }
    return result;
}

std::string MaxGenotypeQuality::do_name() const
{
    return "GQ";
}

std::vector<std::string> MaxGenotypeQuality::do_requirements() const
{
    return {"Samples"};
}

} // namespace csr
} // namespace octopus
