// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "clipped_read_fraction.hpp"

#include <cassert>

#include <boost/variant.hpp>

#include "basics/aligned_read.hpp"
#include "io/variant/vcf_record.hpp"
#include "utils/read_stats.hpp"
#include "../facets/overlapping_reads.hpp"

namespace octopus { namespace csr {

std::unique_ptr<Measure> ClippedReadFraction::do_clone() const
{
    return std::make_unique<ClippedReadFraction>(*this);
}

namespace {

auto clip_fraction(const AlignedRead& read) noexcept
{
    assert(sequence_size(read) > 0);
    return static_cast<double>(total_clip_size(read)) / sequence_size(read);
}

bool is_significantly_clipped(const AlignedRead& read) noexcept
{
    return is_soft_clipped(read) && clip_fraction(read) > 0.25;
}

double clipped_fraction(const ReadMap& reads, const GenomicRegion& region)
{
    unsigned num_reads {0}, num_soft_clipped_reads {0};
    for (const auto& p : reads) {
        for (const auto& read : overlap_range(p.second, region)) {
            if (is_significantly_clipped(read)) ++num_soft_clipped_reads;
            ++num_reads;
        }
    }
    return static_cast<double>(num_soft_clipped_reads) / num_reads;
}

} // namespace

Measure::ResultType ClippedReadFraction::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& reads = get_value<OverlappingReads>(facets.at("OverlappingReads"));
    return clipped_fraction(reads, mapped_region(call));
}

Measure::ResultCardinality ClippedReadFraction::do_cardinality() const noexcept
{
    return ResultCardinality::one;
}

std::string ClippedReadFraction::do_name() const
{
    return "CRF";
}

std::string ClippedReadFraction::do_describe() const
{
    return "Fraction of clipped reads covering the call";
}

std::vector<std::string> ClippedReadFraction::do_requirements() const
{
    return {"OverlappingReads"};
}

} // namespace csr
} // namespace octopus
