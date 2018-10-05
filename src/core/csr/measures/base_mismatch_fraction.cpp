// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "base_mismatch_fraction.hpp"

#include <iterator>
#include <algorithm>
#include <cassert>

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "utils/concat.hpp"
#include "base_mismatch_count.hpp"
#include "variant_length.hpp"

namespace octopus { namespace csr {

const std::string BaseMismatchFraction::name_ = "BMF";

BaseMismatchFraction::BaseMismatchFraction()
: depth_ {true, false}
{}

std::unique_ptr<Measure> BaseMismatchFraction::do_clone() const
{
    return std::make_unique<BaseMismatchFraction>(*this);
}

Measure::ResultType BaseMismatchFraction::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto depths = boost::get<std::vector<std::size_t>>(depth_.evaluate(call, facets));
    const auto mismatch_counts = boost::get<std::vector<int>>(BaseMismatchCount{}.evaluate(call, facets));
    const auto variant_lengths = boost::get<std::vector<int>>(VariantLength{}.evaluate(call, facets));
    assert(depths.size() == mismatch_counts.size());
    assert(mismatch_counts.size() == variant_lengths.size());
    std::vector<int> bases(depths.size());
    std::transform(std::cbegin(depths), std::cend(depths), std::cbegin(variant_lengths), std::begin(bases),
                   [] (const auto& depth, const auto& allele_length) { return depth * allele_length; });
    std::vector<double> result(depths.size());
    std::transform(std::cbegin(mismatch_counts), std::cend(mismatch_counts), std::cbegin(bases), std::begin(result),
                   [] (auto mismatches, auto bases) { return bases > 0 ? static_cast<double>(mismatches) / bases : 0.0; });
    return result;
}

Measure::ResultCardinality BaseMismatchFraction::do_cardinality() const noexcept
{
    return ResultCardinality::num_samples;
}

const std::string& BaseMismatchFraction::do_name() const
{
    return name_;
}

std::string BaseMismatchFraction::do_describe() const
{
    return "Fraction of base mismatches at variant position";
}

std::vector<std::string> BaseMismatchFraction::do_requirements() const
{
    return concat(concat(depth_.requirements(), BaseMismatchCount{}.requirements()), VariantLength{}.requirements());
}
    
} // namespace csr
} // namespace octopus
