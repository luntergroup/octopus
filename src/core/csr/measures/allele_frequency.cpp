// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "allele_frequency.hpp"

#include <algorithm>
#include <iterator>

#include <boost/variant.hpp>

#include "core/types/allele.hpp"
#include "io/variant/vcf_record.hpp"
#include "utils/concat.hpp"
#include "allele_depth.hpp"
#include "assigned_depth.hpp"

namespace octopus { namespace csr {

const std::string AlleleFrequency::name_ = "AF";

std::unique_ptr<Measure> AlleleFrequency::do_clone() const
{
    return std::make_unique<AlleleFrequency>(*this);
}

Measure::ResultType AlleleFrequency::get_default_result() const
{
    return std::vector<boost::optional<double>> {};
}

Measure::ResultType AlleleFrequency::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto allele_depths = boost::get<std::vector<boost::optional<int>>>(AlleleDepth{}.evaluate(call, facets));
    const auto depths = boost::get<std::vector<std::size_t>>(AssignedDepth{}.evaluate(call, facets));
    std::vector<boost::optional<double>> result(allele_depths.size());
    std::transform(std::cbegin(allele_depths), std::cend(allele_depths), std::cbegin(depths), std::begin(result),
                   [] (const auto allele_depth, const auto depth) -> boost::optional<double> {
        if (allele_depth && depth > 0) {
            return static_cast<double>(*allele_depth) / depth;
        } else {
            return boost::none;
        }
    });
    return result;
}

Measure::ResultCardinality AlleleFrequency::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& AlleleFrequency::do_name() const
{
    return name_;
}

std::string AlleleFrequency::do_describe() const
{
    return "Empirical minor allele frequency of ALT alleles (AD / ADP)";
}

std::vector<std::string> AlleleFrequency::do_requirements() const
{
    return concat(AlleleDepth{}.requirements(), AssignedDepth{}.requirements());
}

} // namespace csr
} // namespace octopus
