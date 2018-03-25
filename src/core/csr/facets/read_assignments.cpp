// Copyright (c) 2015-2018 Daniel Cooke
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

bool is_homozygous_nonreference(const Genotype<Haplotype>& genotype)
{
    return genotype.is_homozygous() && !is_reference(genotype[0]);
}

} // namespace

ReadAssignments::ReadAssignments(const ReferenceGenome& reference, const GenotypeMap& genotypes, const ReadMap& reads)
: result_ {}
{
    result_.reserve(genotypes.size());
    for (const auto& p : genotypes) {
        const auto& sample = p.first;
        const auto& genotypes = p.second;
        result_[sample].reserve(genotypes.size());
        for (const auto& genotype : genotypes) {
            auto local_reads = copy_overlapped_to_vector(reads.at(sample), genotype);
            for (const auto& haplotype : genotype) {
                // So every called haplotype appears in support map, even if no read support
                result_[sample][haplotype] = {};
            }
            if (!local_reads.empty()) {
                HaplotypeSupportMap genotype_support {};
                if (!is_homozygous_nonreference(genotype)) {
                    genotype_support = compute_haplotype_support(genotype, local_reads);
                } else {
                    auto augmented_genotype = genotype;
                    Haplotype ref {mapped_region(genotype), reference};
                    result_[sample][ref] = {};
                    augmented_genotype.emplace(std::move(ref));
                    genotype_support = compute_haplotype_support(augmented_genotype, local_reads);
                }
                for (auto& s : genotype_support) {
                    result_[sample][s.first] = std::move(s.second);
                }
            }
        }
    }
}

Facet::ResultType ReadAssignments::do_get() const
{
    return std::cref(result_);
}

} // namespace csr
} // namespace octopus
