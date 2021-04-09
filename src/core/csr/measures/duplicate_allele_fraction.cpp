// Copyright (c) 2015-2021 Daniel Cooke
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

Measure::ValueType DuplicateAlleleFraction::get_value_type() const
{
    return double {};
}

Measure::ResultType DuplicateAlleleFraction::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto allele_depths = boost::get<Array<Array<Optional<ValueType>>>>(AlleleDepth{}.evaluate(call, facets));
    const auto duplicate_allele_depths = boost::get<Array<Array<Optional<ValueType>>>>(DuplicateAlleleDepth{}.evaluate(call, facets));
    const auto num_alleles = call.alt().size() + 1;
    assert(allele_depths.size() == duplicate_allele_depths.size());
    Array<Array<Optional<ValueType>>> result(allele_depths.size(), Array<Optional<ValueType>>(num_alleles));
    for (std::size_t s {0}; s < allele_depths.size(); ++s) {
        assert(allele_depths[s].size() == num_alleles && duplicate_allele_depths[s].size() == num_alleles);
        for (std::size_t a {0}; a < num_alleles; ++a) {
            if (allele_depths[s][a] && duplicate_allele_depths[s][a]) {
                const auto allele_depth = boost::get<std::size_t>(*allele_depths[s][a]);
                if (allele_depth > 0) {
                    result[s][a] = static_cast<double>(boost::get<std::size_t>(*duplicate_allele_depths[s][a])) / allele_depth;
                } else {
                    result[s][a] = 0.0;
                }
            }
        }
    }
    return result;
}

Measure::ResultCardinality DuplicateAlleleFraction::do_cardinality() const noexcept
{
    return ResultCardinality::samples_and_alleles;
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

boost::optional<Measure::Aggregator> DuplicateAlleleFraction::do_aggregator() const noexcept
{
    return Measure::Aggregator::max_tail;
}

} // namespace csr
} // namespace octopus
