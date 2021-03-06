// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "misaligned_read_count.hpp"

#include <algorithm>
#include <iterator>

#include <boost/variant.hpp>

#include "core/tools/read_assigner.hpp"
#include "core/tools/read_realigner.hpp"
#include "core/types/allele.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "../facets/samples.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string MisalignedReadCount::name_ = "MRC";

std::unique_ptr<Measure> MisalignedReadCount::do_clone() const
{
    return std::make_unique<MisalignedReadCount>(*this);
}

Measure::ValueType MisalignedReadCount::get_value_type() const
{
    return int {};
}

double error_expectation(const AlignedRead::BaseQualityVector& qualities)
{
    return std::accumulate(std::cbegin(qualities), std::cend(qualities), 0.0,
                           [] (auto curr, auto q) { return curr + std::pow(10.0, -q / 10); });
}

double error_expectation(const AlignedRead& read)
{
    return error_expectation(read.base_qualities());
}

unsigned count_errors(const AlignedRead& read)
{
    return sum_non_matches(read.cigar());
}

bool is_likely_misaligned(const AlignedRead& read)
{
    const auto expected_errors = static_cast<unsigned>(std::ceil(error_expectation(read)));
    const auto max_expected_errors = 5 * expected_errors;
    const auto observed_errors = count_errors(read);
    return observed_errors > max_expected_errors;
}

std::size_t count_likely_misaligned(const std::vector<AlignedRead>& reads, const GenomicRegion& region)
{
    const auto overlapped = overlap_range(reads, region);
    return std::count_if(std::cbegin(overlapped), std::cend(overlapped), is_likely_misaligned);
}

Measure::ResultType MisalignedReadCount::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).haplotypes;
    Array<ValueType> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        int sample_result {0};
        if (assignments.count(sample) == 1) {
            for (const auto& p : assignments.at(sample).assigned_wrt_haplotype) {
                sample_result += count_likely_misaligned(p.second, mapped_region(call));
            }
        }
        result.emplace_back(sample_result);
    }
    return result;
}

Measure::ResultCardinality MisalignedReadCount::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& MisalignedReadCount::do_name() const
{
    return name_;
}

std::string MisalignedReadCount::do_describe() const
{
    return "Number of reads supporting the call that appear misaligned";
}

std::vector<std::string> MisalignedReadCount::do_requirements() const
{
    return {"Samples", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
