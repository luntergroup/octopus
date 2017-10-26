// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "mapping_quality_divergence.hpp"

#include <cmath>

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "basics/aligned_read.hpp"
#include "../facets/read_assignments.hpp"

#include "utils/read_stats.hpp"

namespace octopus { namespace csr {

std::unique_ptr<Measure> MappingQualityDivergence::do_clone() const
{
    return std::make_unique<MappingQualityDivergence>(*this);
}

Measure::ResultType MappingQualityDivergence::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto assignments = boost::get<ReadAssignments::ResultType>(facets.at("ReadAssignments").get());
    double max_divergence {0};
    for (const auto& p : assignments) {
        if (call.is_heterozygous(p.first)) {
            if (p.second.size() == 2) {
                auto h1_rmq = rmq_mapping_quality(p.second.cbegin()->second);
                auto h2_rmq = rmq_mapping_quality(std::next(p.second.cbegin())->second);
                max_divergence = std::max(max_divergence, std::abs(h1_rmq - h2_rmq));
            }
        }
    }
    return max_divergence;
}

std::string MappingQualityDivergence::do_name() const
{
    return "MQD";
}

std::vector<std::string> MappingQualityDivergence::do_requirements() const
{
    return {"ReadAssignments"};
}

} // namespace csr
} // namespace octopus
