// Copyright (c) 2015-2020 Daniel Cooke
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
    const auto allele_depths = boost::get<Array<Optional<Array<ValueType>>>>(AlleleDepth{}.evaluate(call, facets));
    const auto duplicate_allele_depths = boost::get<Array<Optional<Array<ValueType>>>>(DuplicateAlleleDepth{}.evaluate(call, facets));
    assert(allele_depths.size() == duplicate_allele_depths.size());
    Array<Optional<Array<Optional<ValueType>>>> result(allele_depths.size());
    std::transform(std::cbegin(duplicate_allele_depths), std::cend(duplicate_allele_depths), std::cbegin(allele_depths), std::begin(result),
                   [] (const auto& duplicate_allele_depths, const auto& allele_depths) {
        Optional<Array<Optional<ValueType>>> result {};
        if (duplicate_allele_depths && allele_depths) {
            assert(duplicate_allele_depths->size() == allele_depths->size());
            result = Array<Optional<ValueType>>(allele_depths->size());
            std::transform(std::cbegin(*duplicate_allele_depths), std::cend(*duplicate_allele_depths),
                           std::cbegin(*allele_depths), std::begin(*result),
                           [] (const auto dup_depth, const auto allele_depth) {
                                Optional<ValueType> result {};
                                const auto adp = boost::get<std::size_t>(allele_depth);;
                                if (adp > 0) {
                                    result = static_cast<double>(boost::get<std::size_t>(dup_depth)) / adp;
                                }
                                return result;
                           });
        }
        return result;
    });
    return result;
}

Measure::ResultCardinality DuplicateAlleleFraction::do_cardinality() const noexcept
{
    return ResultCardinality::samples_and_alt_alleles;
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
