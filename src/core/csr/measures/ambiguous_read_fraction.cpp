// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "ambiguous_read_fraction.hpp"

#include <algorithm>
#include <iterator>

#include <boost/variant.hpp>

#include "core/tools/read_assigner.hpp"
#include "core/types/allele.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "utils/genotype_reader.hpp"
#include "../facets/samples.hpp"
#include "../facets/overlapping_reads.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string AmbiguousReadFraction::name_ = "ARF";

std::unique_ptr<Measure> AmbiguousReadFraction::do_clone() const
{
    return std::make_unique<AmbiguousReadFraction>(*this);
}

Measure::ValueType AmbiguousReadFraction::get_value_type() const
{
    return double {};
}

Measure::ResultType AmbiguousReadFraction::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& reads = get_value<OverlappingReads>(facets.at("OverlappingReads"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).haplotypes;
    Array<Optional<ValueType>> result(samples.size());
    for (std::size_t s {0}; s < samples.size(); ++s) {
        const auto& sample = samples[s];
        const auto num_overlapping_reads = count_overlapped(reads.at(sample), call);
        if (num_overlapping_reads > 0) {
            double sample_result {};
            if (assignments.count(sample) == 1) {
                const auto num_ambiguous_reads = count_overlapped(assignments.at(sample).ambiguous_wrt_reference, call);
                sample_result = static_cast<double>(num_ambiguous_reads) / num_overlapping_reads;
            }
            result[s] = sample_result;
        }
    }
    return result;
}

Measure::ResultCardinality AmbiguousReadFraction::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& AmbiguousReadFraction::do_name() const
{
    return name_;
}

std::string AmbiguousReadFraction::do_describe() const
{
    return "Fraction of reads overlapping the call that cannot be assigned to a unique haplotype";
}

std::vector<std::string> AmbiguousReadFraction::do_requirements() const
{
    return {"Samples", "OverlappingReads", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
