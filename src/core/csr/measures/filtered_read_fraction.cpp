// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "filtered_read_fraction.hpp"

#include <boost/variant.hpp>

namespace octopus { namespace csr {

FilteredReadFraction::FilteredReadFraction()
: calling_depth_ {false}
, filtering_depth_ {true}
{}

std::unique_ptr<Measure> FilteredReadFraction::do_clone() const
{
    return std::make_unique<FilteredReadFraction>(*this);
}

Measure::ResultType FilteredReadFraction::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto filtering_depth = boost::get<std::size_t>(filtering_depth_.evaluate(call, facets));
    double result {0};
    if (filtering_depth > 0) {
        auto calling_depth = boost::get<std::size_t>(calling_depth_.evaluate(call, facets));
        result = 1.0 - (static_cast<double>(calling_depth) / filtering_depth);
    }
    return result;
}

Measure::ResultCardinality FilteredReadFraction::do_cardinality() const noexcept
{
    return ResultCardinality::one;
}

std::string FilteredReadFraction::do_name() const
{
    return "FRF";
}

std::vector<std::string> FilteredReadFraction::do_requirements() const
{
    return filtering_depth_.requirements();
}

} // namespace csr
} // namespace octopus
