// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "strand_disequilibrium.hpp"

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

namespace {

bool includes(const std::pair<double, double>& credible_region, const double x)
{
    return credible_region.first <= x && x <= credible_region.second;
}

auto calculate_required_balance_mass(const double a, const double b, const double tolerance = 0.01)
{
    // TODO: can we do this analytically?
    double result {tolerance};
    while (result <= 1.0) {
        const auto credible_region = maths::beta_hdi(a, b, result);
        if (includes(credible_region, 0.5)) break;
        result *= 2;
    }
    return std::min(result, 1.0);
}

auto calculate_distance_to_half(const std::pair<double, double>& credible_region)
{
    if (includes(credible_region, 0.5)) {
        return 0.0;
    } else if (credible_region.second < 0.5) {
        return 0.5 - credible_region.second;
    } else {
        return credible_region.second - 0.5;
    }
}

} // namespace

Measure::ResultType StrandDisequilibrium::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& reads = get_value<OverlappingReads>(facets.at("OverlappingReads"));
    std::vector<double> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        const auto direction_counts = count_directions(reads.at(sample), mapped_region(call));
        const auto credible_region = maths::beta_hdi(direction_counts.first + 0.5, direction_counts.second + 0.5);
        result.push_back(calculate_distance_to_half(credible_region));
    }
    return result;
}

Measure::ResultCardinality StrandDisequilibrium::do_cardinality() const noexcept
{
    return ResultCardinality::num_samples;
}

const std::string& StrandDisequilibrium::do_name() const
{
    return name_;
}

std::string StrandDisequilibrium::do_describe() const
{
    return "Strand bias of reads overlapping the site; probability mass required to include 0.5 in credible interval";
}

std::vector<std::string> StrandDisequilibrium::do_requirements() const
{
    return {"Samples", "OverlappingReads"};
}

} // namespace csr
} // namespace octopus