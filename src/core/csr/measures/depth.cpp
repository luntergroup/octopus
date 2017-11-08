// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "depth.hpp"

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "../facets/overlapping_reads.hpp"

namespace octopus { namespace csr {

Depth::Depth(bool recalculate) : recalculate_ {recalculate} {}

std::unique_ptr<Measure> Depth::do_clone() const
{
    return std::make_unique<Depth>(*this);
}

Measure::ResultType Depth::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    if (recalculate_) {
        const auto reads = boost::get<OverlappingReads::ResultType>(facets.at("OverlappingReads").get());
        return static_cast<std::size_t>(count_overlapped(reads, call));
    } else {
        return static_cast<std::size_t>(std::stoull(call.info_value(vcfspec::info::combinedReadDepth).front()));
    }
}

std::string Depth::do_name() const
{
    return "DP";
}

std::vector<std::string> Depth::do_requirements() const
{
    if (recalculate_) {
        return {"OverlappingReads"};
    } else {
        return {};
    }
}

} // namespace csr
} // namespace octopus