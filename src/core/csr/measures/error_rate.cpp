// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "error_rate.hpp"

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
#include "../facets/samples.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string ErrorRate::name_ = "ER";

std::unique_ptr<Measure> ErrorRate::do_clone() const
{
    return std::make_unique<ErrorRate>(*this);
}

namespace {

boost::optional<double> 
compute_error_rate(const Facet::SupportMaps& assignments, const SampleName& sample) noexcept
{
    std::size_t error_bases {0}, total_bases {0};
    if (assignments.support.count(sample) == 1) {
        for (const auto& p : assignments.support.at(sample)) {
            for (auto read : safe_realign(p.second, p.first)) {
                error_bases += sum_non_matches(read.cigar());
                total_bases += sequence_size(read);
            }
        }
    }
    if (assignments.ambiguous.count(sample) == 1) {
        std::map<Haplotype, std::vector<AlignedRead>> assigned {};
        for (const auto& ambiguous_read : assignments.ambiguous.at(sample)) {
            if (ambiguous_read.haplotypes) {
                assigned[*ambiguous_read.haplotypes->front()].push_back(ambiguous_read.read);
            }
        }
        for (auto& p : assigned) {
            safe_realign(p.second, p.first);
            for (auto read : p.second) {
                error_bases += sum_non_matches(read.cigar());
                total_bases += sequence_size(read);
            }
        }
    }
    if (total_bases > 0) {
        return static_cast<double>(error_bases) / total_bases;
    } else {
        return boost::none;
    }
}

} // namespace

Measure::ResultType ErrorRate::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments"));
    std::vector<boost::optional<double>> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        result.push_back(compute_error_rate(assignments, sample));
    }
    return result;
}

Measure::ResultCardinality ErrorRate::do_cardinality() const noexcept
{
    return ResultCardinality::num_samples;
}

const std::string& ErrorRate::do_name() const
{
    return name_;
}

std::string ErrorRate::do_describe() const
{
    return "Error rate in reads overlapping the site";
}

std::vector<std::string> ErrorRate::do_requirements() const
{
    return {"Samples", "ReadAssignments"};
}
    
} // namespace csr
} // namespace octopus
