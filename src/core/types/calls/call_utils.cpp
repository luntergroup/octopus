// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "call_utils.hpp"

namespace octopus {

namespace detail {

std::vector<unsigned>
compute_haplotype_order(const std::vector<Genotype<Allele>>& genotypes, const ReferenceGenome& reference)
{
    if (!genotypes.empty()) {
        const auto ploidy = genotypes.front().ploidy();
        const auto region = encompassing_region(genotypes);
        assert(std::all_of(std::cbegin(genotypes), std::cend(genotypes), [ploidy] (const auto& g) { return g.ploidy() == ploidy; }));
        std::vector<std::pair<Haplotype, unsigned>> haplotypes {};
        haplotypes.reserve(ploidy);
        for (unsigned haplotype_idx {0}; haplotype_idx < ploidy; ++haplotype_idx) {
            Haplotype::Builder haplotype {region, reference};
            for (const auto& genotype : genotypes) {
                const auto& allele = genotype[haplotype_idx];
                if (haplotype.can_push_back(allele)) {
                    haplotype.push_back(allele);
                }
            }
            haplotypes.push_back({haplotype.build(), haplotype_idx});
        }
        std::sort(std::begin(haplotypes), std::end(haplotypes));
        std::vector<unsigned> result(ploidy);
        std::transform(std::cbegin(haplotypes), std::cend(haplotypes), std::begin(result),
                       [] (const auto& p) { return p.second; });
        return result;
    } else {
        return {};
    }
}

bool are_same_phase_set(const Call::GenotypeCall& lhs, const Call::GenotypeCall& rhs)
{
    return lhs.phase && rhs.phase && overlaps(lhs.phase->region(), rhs.phase->region());
}

bool are_same_phase_set(const CallWrapper& lhs, const CallWrapper& rhs, const SampleName& sample)
{
    return are_same_phase_set(lhs->get_genotype_call(sample), rhs->get_genotype_call(sample));
}

bool have_different_phase_sets(const Call::GenotypeCall& lhs, const Call::GenotypeCall& rhs)
{
    return lhs.phase && rhs.phase && !overlaps(lhs.phase->region(), rhs.phase->region());
}

bool have_different_phase_sets(const CallWrapper& lhs, const CallWrapper& rhs, const SampleName& sample)
{
    return have_different_phase_sets(lhs->get_genotype_call(sample), rhs->get_genotype_call(sample));
}

} // namespace detail

} // namespace octopus

