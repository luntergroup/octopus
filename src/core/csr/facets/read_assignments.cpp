// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_assignments.hpp"

namespace octopus { namespace csr {

const std::string ReadAssignments::name_ {"ReadAssignments"};

namespace {

template <typename Mappable>
auto copy_overlapped_to_vector(const ReadContainer& reads, const Mappable& mappable)
{
    const auto overlapped = overlap_range(reads, mappable);
    return std::vector<AlignedRead> {std::cbegin(overlapped), std::cend(overlapped)};
}

} // namespace

ReadAssignments::ReadAssignments(const GenotypeMap& genotypes, const ReadMap& reads)
: assignments_ {}
{
    assignments_.reserve(genotypes.size());
    for (const auto& p : genotypes) {
        assignments_[p.first].reserve(p.second.size());
        for (const auto& genotype : p.second) {
            auto local_reads = copy_overlapped_to_vector(reads.at(p.first), genotype);
            for (const auto& haplotype : genotype) {
                // So every called haplotype appears in support map, even if no read support
                assignments_[p.first][haplotype] = {};
            }
            if (!local_reads.empty()) {
                auto local_assignments = compute_haplotype_support(genotype, local_reads);
                for (auto s : local_assignments) {
                    assignments_[p.first][s.first] = s.second;
                }
            }
        }
    }
}

Facet::ResultType ReadAssignments::do_get() const
{
    return assignments_;
}

} // namespace csr
} // namespace octopus
