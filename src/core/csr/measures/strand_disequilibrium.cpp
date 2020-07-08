// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "strand_disequilibrium.hpp"

#include <boost/lexical_cast.hpp>

#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "basics/aligned_read.hpp"
#include "utils/read_stats.hpp"
#include "utils/maths.hpp"
#include "../facets/samples.hpp"
#include "../facets/overlapping_reads.hpp"

namespace octopus { namespace csr {

const std::string StrandDisequilibrium::name_ = "SD";

std::unique_ptr<Measure> StrandDisequilibrium::do_clone() const
{
    return std::make_unique<StrandDisequilibrium>(*this);
}

Measure::ValueType StrandDisequilibrium::get_value_type() const
{
    return double {};
}

void StrandDisequilibrium::do_set_parameters(std::vector<std::string> params)
{
    if (params.size() != 1) {
        throw BadMeasureParameters {this->name(), "only has one parameter (tail mass)"};
    }
    try {
        tail_mass_ = boost::lexical_cast<decltype(tail_mass_)>(params.front());
    } catch (const boost::bad_lexical_cast&) {
        throw BadMeasureParameters {this->name(), "given parameter \"" + params.front() + "\" cannot be parsed"};
    }
    if (tail_mass_ < 0 || tail_mass_ > 1) {
        throw BadMeasureParameters {this->name(), "tail mass must be between 0 and 1"};
    }
}

std::vector<std::string> StrandDisequilibrium::do_parameters() const
{
    return {utils::to_string(tail_mass_, 2)};
}

Measure::ResultType StrandDisequilibrium::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& reads = get_value<OverlappingReads>(facets.at("OverlappingReads"));
    Array<ValueType> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        const auto direction_counts = count_directions(reads.at(sample), mapped_region(call));
        const auto tail_probability = maths::beta_tail_probability(direction_counts.first + 0.5, direction_counts.second + 0.5, tail_mass_);
        result.emplace_back(tail_probability);
    }
    return result;
}

Measure::ResultCardinality StrandDisequilibrium::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& StrandDisequilibrium::do_name() const
{
    return name_;
}

std::string StrandDisequilibrium::do_describe() const
{
    return "Strand bias of reads overlapping the site; probability mass in tails of Beta distribution";
}

std::vector<std::string> StrandDisequilibrium::do_requirements() const
{
    return {"Samples", "OverlappingReads"};
}

bool StrandDisequilibrium::is_equal(const Measure& other) const noexcept
{
    return tail_mass_ == static_cast<const StrandDisequilibrium&>(other).tail_mass_;
}

} // namespace csr
} // namespace octopus
