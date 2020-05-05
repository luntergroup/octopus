// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "is_transversion.hpp"

#include <iterator>
#include <algorithm>

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "core/types/variant.hpp"

namespace octopus { namespace csr {

const std::string IsTransversion::name_ = "ITV";

std::unique_ptr<Measure> IsTransversion::do_clone() const
{
    return std::make_unique<IsTransversion>(*this);
}

Measure::ValueType IsTransversion::get_value_type() const
{
    return bool {};
}

Measure::ResultType IsTransversion::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    Array<ValueType> result(call.alt().size());
    const auto is_transversion_helper = [&] (const auto& alt) {
        const Variant variant {call.mapped_region(), call.ref(), alt};
        return is_transversion(variant);
    };
    std::transform(std::cbegin(call.alt()), std::cend(call.alt()), std::begin(result), is_transversion_helper);
    return result;
}

Measure::ResultCardinality IsTransversion::do_cardinality() const noexcept
{
    return ResultCardinality::alt_alleles;
}

const std::string& IsTransversion::do_name() const
{
    return name_;
}

std::string IsTransversion::do_describe() const
{
    return "Is the variant a transversion";
}

std::vector<std::string> IsTransversion::do_requirements() const
{
    return {};
}

} // namespace csr
} // namespace octopus
