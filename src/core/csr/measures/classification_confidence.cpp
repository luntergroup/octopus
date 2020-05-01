// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "classification_confidence.hpp"

#include "io/variant/vcf_record.hpp"

#include "quality.hpp"
#include "posterior_probability.hpp"

namespace octopus { namespace csr {

const std::string ClassificationConfidence::name_ = "CC";

std::unique_ptr<Measure> ClassificationConfidence::do_clone() const
{
    return std::make_unique<ClassificationConfidence>(*this);
}

Measure::ResultType ClassificationConfidence::get_default_result() const
{
    return boost::optional<double> {};
}

Measure::ResultType ClassificationConfidence::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto quality = boost::get<boost::optional<double>>(Quality().evaluate(call, facets));
    const auto posterior = boost::get<boost::optional<double>>(PosteriorProbability().evaluate(call, facets));
    boost::optional<double> result {};
    if (quality && posterior) {
        if (*quality > 0) {
            result = *posterior / *quality;
        } else {
            result = 0;
        }
    }
    return result;
}

Measure::ResultCardinality ClassificationConfidence::do_cardinality() const noexcept
{
    return ResultCardinality::one;
}

const std::string& ClassificationConfidence::do_name() const
{
    return name_;
}

std::string ClassificationConfidence::do_describe() const
{
    return "PP divided by QUAL";
}

} // namespace csr
} // namespace octopus
