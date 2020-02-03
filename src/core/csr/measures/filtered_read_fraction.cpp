// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "filtered_read_fraction.hpp"

#include <iterator>
#include <algorithm>

#include <boost/variant.hpp>

namespace octopus { namespace csr {

const std::string FilteredReadFraction::name_ = "FRF";

FilteredReadFraction::FilteredReadFraction(bool aggregate_samples)
: calling_depth_ {false, aggregate_samples}
, filtering_depth_ {true, aggregate_samples}
{}

std::unique_ptr<Measure> FilteredReadFraction::do_clone() const
{
    return std::make_unique<FilteredReadFraction>(*this);
}

Measure::ResultType FilteredReadFraction::get_default_result() const
{
    return std::vector<double> {};
}

Measure::ResultType FilteredReadFraction::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    if (filtering_depth_.cardinality() == Measure::ResultCardinality::samples) {
        const auto calling_depth   = boost::get<std::vector<std::size_t>>(calling_depth_.evaluate(call, facets));
        const auto filtering_depth = boost::get<std::vector<std::size_t>>(filtering_depth_.evaluate(call, facets));
        assert(calling_depth.size() == filtering_depth.size());
        std::vector<double> result(calling_depth.size());
        std::transform(std::cbegin(calling_depth), std::cend(calling_depth), std::cbegin(filtering_depth), std::begin(result),
                       [] (auto cd, auto fd) { return fd > 0 ? 1.0 - (static_cast<double>(cd) / fd) : 0; });
        return result;
    } else {
        const auto filtering_depth = boost::get<std::size_t>(filtering_depth_.evaluate(call, facets));
        double result {0};
        if (filtering_depth > 0) {
            auto calling_depth = boost::get<std::size_t>(calling_depth_.evaluate(call, facets));
            result = 1.0 - (static_cast<double>(calling_depth) / filtering_depth);
        }
        return result;
    }
}

Measure::ResultCardinality FilteredReadFraction::do_cardinality() const noexcept
{
    return filtering_depth_.cardinality();
}

const std::string& FilteredReadFraction::do_name() const
{
    return name_;
}

std::string FilteredReadFraction::do_describe() const
{
    return "Fraction of reads filtered for calling";
}

std::vector<std::string> FilteredReadFraction::do_requirements() const
{
    return filtering_depth_.requirements();
}

bool FilteredReadFraction::is_equal(const Measure& other) const noexcept
{
    return calling_depth_ == static_cast<const FilteredReadFraction&>(other).calling_depth_;
}

} // namespace csr
} // namespace octopus
