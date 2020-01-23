// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "duplicate_allele_fraction.hpp"

#include <algorithm>
#include <iterator>
#include <cassert>

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "allele_depth.hpp"
#include "duplicate_allele_depth.hpp"

namespace octopus { namespace csr {

const std::string DuplicateAlleleFraction::name_ = "DAF";

std::unique_ptr<Measure> DuplicateAlleleFraction::do_clone() const
{
    return std::make_unique<DuplicateAlleleFraction>(*this);
}

Measure::ResultType DuplicateAlleleFraction::get_default_result() const
{
    return std::vector<boost::optional<double>> {};
}

Measure::ResultType DuplicateAlleleFraction::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto allele_depths = boost::get<std::vector<boost::optional<int>>>(AlleleDepth{}.evaluate(call, facets));
    const auto duplicate_allele_depths = boost::get<std::vector<boost::optional<int>>>(DuplicateAlleleDepth{}.evaluate(call, facets));
    assert(allele_depths.size() == duplicate_allele_depths.size());
    std::vector<boost::optional<double>> result(allele_depths.size());
    std::transform(std::cbegin(allele_depths), std::cend(allele_depths), std::cbegin(duplicate_allele_depths), std::begin(result),
                   [] (auto depth, auto duplicate_depth) -> boost::optional<double> {
        if (depth && duplicate_depth) {
            return *depth > 0 ? static_cast<double>(std::min(*duplicate_depth, *depth)) / *depth : 0;
        } else {
            return boost::none;
        }
    });
    return result;
}

Measure::ResultCardinality DuplicateAlleleFraction::do_cardinality() const noexcept
{
    return ResultCardinality::num_samples;
}

const std::string& DuplicateAlleleFraction::do_name() const
{
    return name_;
}

std::string DuplicateAlleleFraction::do_describe() const
{
    return "Fraction of realigned reads supporting ALT alleles identified as duplicates";
}

std::vector<std::string> DuplicateAlleleFraction::do_requirements() const
{
    
    return {"Samples", "OverlappingReads", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
