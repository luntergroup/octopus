// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "is_somatic.hpp"

#include <iterator>
#include <algorithm>
#include <cassert>

#include "io/variant/vcf_record.hpp"
#include "../facets/samples.hpp"
#include "../facets/ploidies.hpp"

namespace octopus { namespace csr {

const std::string IsSomatic::name_ = "SOMATIC";

IsSomatic::IsSomatic(bool report_sample_status) : report_sample_status_ {report_sample_status} {}

std::unique_ptr<Measure> IsSomatic::do_clone() const
{
    return std::make_unique<IsSomatic>(*this);
}

Measure::ValueType IsSomatic::get_value_type() const
{
    return bool {};
}

namespace {

bool is_somatic_sample(const VcfRecord& call, const VcfRecord::SampleName& sample, const unsigned sample_ploidy)
{
    assert(is_somatic(call));
    const auto observed_sample_ploidy = call.ploidy(sample);
    if (observed_sample_ploidy > sample_ploidy) {
        return call.has_alt_allele(sample);
    } else {
        return false;
    }
}

} // namespace

Measure::ResultType IsSomatic::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    if (report_sample_status_) {
        const auto& samples = get_value<Samples>(facets.at("Samples"));
        Array<ValueType> result(samples.size(), false);
        if (is_somatic(call)) {
            const auto& ploidies = get_value<Ploidies>(facets.at("Ploidies"));
            std::transform(std::cbegin(samples), std::cend(samples), std::begin(result),
                           [&] (const auto& sample) { return is_somatic_sample(call, sample, ploidies.at(sample)); });
        }
        return result;
    } else {
        return ValueType {is_somatic(call)};
    }
}

Measure::ResultCardinality IsSomatic::do_cardinality() const noexcept
{
    if (report_sample_status_) {
        return ResultCardinality::samples;
    } else {
        return ResultCardinality::one;
    }
}

const std::string& IsSomatic::do_name() const
{
    return name_;
}

std::string IsSomatic::do_describe() const
{
    if (report_sample_status_) {
        return "SOMATIC status of each sample";
    } else {
        return "Is the call marked SOMATIC";
    }
}

std::vector<std::string> IsSomatic::do_requirements() const
{
    if (report_sample_status_) {
        return {"Samples", "Ploidies"};
    } else {
        return {};
    }
}

bool IsSomatic::is_equal(const Measure& other) const noexcept
{
    return report_sample_status_ == static_cast<const IsSomatic&>(other).report_sample_status_;
}

} // namespace csr
} // namespace octopus
