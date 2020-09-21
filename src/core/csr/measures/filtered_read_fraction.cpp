// Copyright (c) 2015-2020 Daniel Cooke
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

Measure::ValueType FilteredReadFraction::get_value_type() const
{
    return double {};
}

Measure::ResultType FilteredReadFraction::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    if (filtering_depth_.cardinality() == Measure::ResultCardinality::samples) {
        const auto raw_depth = boost::get<Array<ValueType>>(filtering_depth_.evaluate(call, facets));
        const auto filtered_depth   = boost::get<Array<ValueType>>(calling_depth_.evaluate(call, facets));
        assert(raw_depth.size() == filtered_depth.size());
        Array<ValueType> result(raw_depth.size());
        const auto filtered_read_fraction_helper = [] (const auto raw_depth, const auto filtered_depth) {
            return raw_depth > 0 ? 1.0 - (static_cast<double>(filtered_depth) / raw_depth) : 0.0;
        };
        const auto filtered_read_fraction = [&] (const auto raw_depth, const auto filtered_depth) -> ValueType {
            return filtered_read_fraction_helper(boost::get<std::size_t>(raw_depth), boost::get<std::size_t>(filtered_depth));
        };
        std::transform(std::cbegin(raw_depth), std::cend(raw_depth), std::cbegin(filtered_depth),
                       std::begin(result), filtered_read_fraction);
        return result;
    } else {
        const auto raw_depth = boost::get<std::size_t>(boost::get<ValueType>(filtering_depth_.evaluate(call, facets)));
        ValueType result {0};
        if (raw_depth > 0) {
            auto filtered_depth = boost::get<std::size_t>(boost::get<ValueType>(calling_depth_.evaluate(call, facets)));
            result = 1.0 - (static_cast<double>(filtered_depth) / raw_depth);
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
