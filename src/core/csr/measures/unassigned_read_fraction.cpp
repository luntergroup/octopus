// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "unassigned_read_fraction.hpp"

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

std::unique_ptr<Measure> UnassignedReadFraction::do_clone() const
{
    return std::make_unique<UnassignedReadFraction>(*this);
}

Measure::ResultType UnassignedReadFraction::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& reads = get_value<OverlappingReads>(facets.at("OverlappingReads"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments"));
    boost::optional<double> result {};
    for (const auto& sample : samples) {
        const auto num_overlapping_reads = count_overlapped(reads.at(sample), call);
        if (num_overlapping_reads > 0) {
            std::vector<Allele> alleles; bool has_ref;
            std::tie(alleles, has_ref) = get_called_alleles(call, sample);
            std::vector<AlignedRead> assigned_reads {};
            const auto& support = assignments.at(sample);
            for (const auto& allele : alleles) {
                std::vector<AlignedRead> allele_realignments {};
                for (const auto& h : support) {
                    const auto& haplotype = h.first;
                    const auto& supporting_reads = h.second;
                    if (!supporting_reads.empty() && haplotype.includes(allele)) {
                        for (const auto& read : supporting_reads) {
                            if (overlaps(read, allele)) {
                                assigned_reads.push_back(read);
                            }
                        }
                    }
                }
            }
            const auto assigned_fraction = static_cast<double>(assigned_reads.size()) / num_overlapping_reads;
            if (result) {
                result = std::max(*result, 1 - assigned_fraction);
            } else {
                result = 1 - assigned_fraction;
            }
        }
    }
    return result;
}

Measure::ResultCardinality UnassignedReadFraction::do_cardinality() const noexcept
{
    return ResultCardinality::one;
}

std::string UnassignedReadFraction::do_name() const
{
    return "URF";
}

std::string UnassignedReadFraction::do_describe() const
{
    return "Fraction of reads overlapping the call that cannot be assigned to a haplotype";
}

std::vector<std::string> UnassignedReadFraction::do_requirements() const
{
    return {"Samples", "OverlappingReads", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
