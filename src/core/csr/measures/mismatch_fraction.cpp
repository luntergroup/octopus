// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "mismatch_fraction.hpp"

#include <iterator>
#include <algorithm>
#include <cassert>

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "utils/concat.hpp"

namespace octopus { namespace csr {

const std::string MismatchFraction::name_ = "MF";

MismatchFraction::MismatchFraction()
: mismatch_count_ {}
, depth_ {true, false}
{}

std::unique_ptr<Measure> MismatchFraction::do_clone() const
{
    return std::make_unique<MismatchFraction>(*this);
}

Measure::ValueType MismatchFraction::get_value_type() const
{
    return double {};
}

Measure::ResultType MismatchFraction::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto depths = boost::get<Array<ValueType>>(depth_.evaluate(call, facets));
    const auto mismatch_counts = boost::get<Array<ValueType>>(mismatch_count_.evaluate(call, facets));
    assert(depths.size() == mismatch_counts.size());
    Array<ValueType> result(depths.size());
    const auto mismatch_fraction_helper = [] (const auto mismatches, const auto depth) {
        return depth > 0 ? static_cast<double>(mismatches) / depth : 0.0;
    };
    const auto mismatch_fraction = [&] (const auto mismatches, const auto depth) {
        return mismatch_fraction_helper(boost::get<int>(mismatches), boost::get<std::size_t>(depth));
    };
    std::transform(std::cbegin(mismatch_counts), std::cend(mismatch_counts), std::cbegin(depths),
                   std::begin(result), mismatch_fraction);
    return result;
}

Measure::ResultCardinality MismatchFraction::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& MismatchFraction::do_name() const
{
    return name_;
}

std::string MismatchFraction::do_describe() const
{
    return "Fraction of reads with mismatches at variant position";
}

std::vector<std::string> MismatchFraction::do_requirements() const
{
    return concat(depth_.requirements(), mismatch_count_.requirements());
}

} // namespace csr
} // namespace octopus
