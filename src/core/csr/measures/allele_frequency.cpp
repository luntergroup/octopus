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

Measure::ValueType AlleleFrequency::get_value_type() const
{
    return double {};
}

Measure::ResultType AlleleFrequency::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto allele_depths = boost::get<Array<Optional<Array<ValueType>>>>(AlleleDepth{}.evaluate(call, facets));
    const auto depths = boost::get<Array<ValueType>>(AssignedDepth{}.evaluate(call, facets));
    Array<Optional<Array<ValueType>>> result(allele_depths.size());
    for (std::size_t s {0}; s < allele_depths.size(); ++s) {
        const auto depth = boost::get<std::size_t>(depths[s]);
        if (allele_depths[s] && depth > 0) {
            result[s] = Array<ValueType>(allele_depths[s]->size());
            std::transform(std::cbegin(*allele_depths[s]), std::cend(*allele_depths[s]), std::begin(*result[s]),
                           [=] (auto allele_depth) {
                               return static_cast<double>(boost::get<std::size_t>(allele_depth)) / depth;
                           });
        }
    }
    return result;
}

Measure::ResultCardinality AlleleFrequency::do_cardinality() const noexcept
{
    return ResultCardinality::samples_and_alt_alleles;
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
