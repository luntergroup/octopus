// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "error_rate_stdev.hpp"

#include <numeric>
#include <algorithm>
#include <iterator>
#include <map>
#include <cstddef>

#include <boost/variant.hpp>

#include "basics/cigar_string.hpp"
#include "basics/aligned_read.hpp"
#include "core/types/haplotype.hpp"
#include "io/variant/vcf_record.hpp"
#include "core/tools/read_realigner.hpp"
#include "utils/maths.hpp"
#include "../facets/samples.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string ErrorRateStdev::name_ = "ERS";

std::unique_ptr<Measure> ErrorRateStdev::do_clone() const
{
    return std::make_unique<ErrorRateStdev>(*this);
}

namespace {

auto calculate_error_rate(const AlignedRead& read) noexcept
{
    return static_cast<double>(sum_non_matches(read.cigar())) / sequence_size(read);
}

boost::optional<double>
compute_error_rate_stdev(const Facet::SupportMaps& assignments, const SampleName& sample, const GenomicRegion& region) noexcept
{
    std::vector<double> error_rates {};
    if (assignments.support.count(sample) == 1) {
        for (const auto& p : assignments.support.at(sample)) {
            auto realigned_reads = copy_overlapped(p.second, region);
            safe_realign(realigned_reads, p.first);
            error_rates.reserve(error_rates.size() + realigned_reads.size());
            for (const auto& read : realigned_reads) {
                error_rates.push_back(calculate_error_rate(read));
            }
        }
    }
    if (assignments.ambiguous.count(sample) == 1) {
        std::map<Haplotype, std::vector<AlignedRead>> assigned {};
        std::size_t num_ambiguous_reads {0};
        for (const auto& ambiguous_read : assignments.ambiguous.at(sample)) {
            if (ambiguous_read.haplotypes && overlaps(ambiguous_read.read, region)) {
                assigned[*ambiguous_read.haplotypes->front()].push_back(ambiguous_read.read);
                ++num_ambiguous_reads;
            }
        }
        error_rates.reserve(error_rates.size() + num_ambiguous_reads);
        for (auto& p : assigned) {
            safe_realign(p.second, p.first);
            for (const auto& read : p.second) {
                error_rates.push_back(calculate_error_rate(read));
            }
        }
    }
    if (!error_rates.empty()) {
        return maths::stdev(error_rates);
    } else {
        return boost::none;
    }
}

} // namespace

Measure::ResultType ErrorRateStdev::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments"));
    std::vector<boost::optional<double>> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        result.push_back(compute_error_rate_stdev(assignments, sample, mapped_region(call)));
    }
    return result;
}

Measure::ResultCardinality ErrorRateStdev::do_cardinality() const noexcept
{
    return ResultCardinality::num_samples;
}

const std::string& ErrorRateStdev::do_name() const
{
    return name_;
}

std::string ErrorRateStdev::do_describe() const
{
    return "Error rate in reads overlapping the site";
}

std::vector<std::string> ErrorRateStdev::do_requirements() const
{
    return {"Samples", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
