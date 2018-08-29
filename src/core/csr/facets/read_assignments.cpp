// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_assignments.hpp"

#include "core/tools/read_realigner.hpp"

namespace octopus { namespace csr {

const std::string ReadAssignments::name_ {"ReadAssignments"};

namespace {

template <typename Mappable>
auto copy_overlapped_to_vector(const ReadContainer& reads, const Mappable& mappable)
{
    const auto overlapped = overlap_range(reads, mappable);
    return std::vector<AlignedRead> {std::cbegin(overlapped), std::cend(overlapped)};
}

void move_insert(std::deque<AlignedRead>& reads, const SampleName& sample, ReadMap& result)
{
    result[sample].insert(std::make_move_iterator(std::begin(reads)), std::make_move_iterator(std::end(reads)));
}

} // namespace

ReadAssignments::ReadAssignments(const ReferenceGenome& reference, const GenotypeMap& genotypes, const ReadMap& reads)
: result_ {}
{
    const auto num_samples = genotypes.size();
    result_.support.reserve(num_samples);
    result_.ambiguous.reserve(num_samples);
    for (const auto& p : genotypes) {
        const auto& sample = p.first;
        const auto& sample_genotypes = p.second;
        result_.support[sample].reserve(sample_genotypes.size());
        for (const auto& genotype : sample_genotypes) {
            auto local_reads = copy_overlapped_to_vector(reads.at(sample), genotype);
            for (const auto& haplotype : genotype) {
                // So every called haplotype appears in support map, even if no read support
                result_.support[sample][haplotype] = {};
            }
            if (!local_reads.empty()) {
                HaplotypeSupportMap genotype_support {};
                std::deque<AlignedRead> unassigned {};
                if (!genotype.is_homozygous()) {
                    genotype_support = compute_haplotype_support(genotype, local_reads, unassigned);
                } else {
                    if (is_reference(genotype[0])) {
                        genotype_support[genotype[0]] = std::move(local_reads);
                    } else {
                        auto augmented_genotype = genotype;
                        Haplotype ref {mapped_region(genotype), reference};
                        result_.support[sample][ref] = {};
                        augmented_genotype.emplace(std::move(ref));
                        genotype_support = compute_haplotype_support(augmented_genotype, local_reads, unassigned);
                    }
                }
                for (auto& s : genotype_support) {
                    safe_realign_to_reference(s.second, s.first);
                    result_.support[sample][s.first] = std::move(s.second);
                }
                move_insert(unassigned, sample, result_.ambiguous);
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
