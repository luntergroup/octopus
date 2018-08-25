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

Measure::ResultType StrandDisequilibrium::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& reads = get_value<OverlappingReads>(facets.at("OverlappingReads"));
    std::vector<double> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        const auto direction_counts = count_directions(reads.at(sample), mapped_region(call));
        const auto tail_probability = maths::beta_tail_probability(direction_counts.first + 0.5, direction_counts.second + 0.5, tail_mass_);
        result.push_back(tail_probability);
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
    return "Strand bias of reads overlapping the site; probability mass in tails of Beta distribution";
}

std::vector<std::string> StrandDisequilibrium::do_requirements() const
{
    return {"Samples", "OverlappingReads"};
}

} // namespace csr
} // namespace octopus