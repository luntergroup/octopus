// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef call_utils_hpp
#define call_utils_hpp

#include <vector>
#include <iterator>
#include <cstddef>
#include <algorithm>

#include "core/types/genotype.hpp"
#include "core/types/haplotype.hpp"
#include "io/reference/reference_genome.hpp"
#include "config/common.hpp"
#include "call_wrapper.hpp"

namespace octopus {

namespace detail {

struct AlleleExtractor
{
    explicit AlleleExtractor(unsigned index) noexcept : index_ {index} {}
    const Allele& operator()(const Genotype<Allele>& genotype) const noexcept { return genotype[index_]; }
private:
    unsigned index_;
};

void sort_alleles_in_haplotype_order(std::vector<Genotype<Allele>>& genotypes, const ReferenceGenome& reference)
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
                if (!is_reference(allele, reference) && haplotype.can_push_back(allele)) {
                    haplotype.push_back(allele);
                }
            }
            haplotypes.push_back({haplotype.build(), haplotype_idx});
        }
        std::sort(std::begin(haplotypes), std::end(haplotypes));
        std::vector<unsigned> allele_order(ploidy);
        std::transform(std::cbegin(haplotypes), std::cend(haplotypes), std::begin(allele_order),
                       [] (const auto& p) { return p.second; });
        if (!std::is_sorted(std::cbegin(allele_order), std::cend(allele_order))) {
            for (auto& genotype : genotypes) {
                genotype.reorder_alleles(allele_order);
            }
        }
    }
}

} // namespace detail

template <typename Iterator>
void sort_genotype_alleles(Iterator first_call, Iterator last_call, const SampleName& sample, const ReferenceGenome& reference)
{
    const auto num_calls = static_cast<std::size_t>(std::distance(first_call, last_call));
    if (num_calls > 0) {
        std::vector<Genotype<Allele>> genotypes {};
        genotypes.reserve(num_calls);
        std::for_each(first_call, last_call, [&] (CallWrapper& call) {
            genotypes.push_back(call->get_genotype_call(sample).genotype);
        });
        detail::sort_alleles_in_haplotype_order(genotypes, reference);
        std::size_t g {0};
        std::for_each(first_call, last_call, [&] (CallWrapper& call) {
            call->get_genotype_call(sample).genotype = std::move(genotypes[g++]);
        });
    }
}

namespace detail {

bool are_in_phase(const Call::GenotypeCall& lhs, const Call::GenotypeCall& rhs)
{
    return lhs.phase && rhs.phase && overlaps(lhs.phase->region(), rhs.phase->region());
}

bool are_in_phase(const CallWrapper& lhs, const CallWrapper& rhs, const SampleName& sample)
{
    return are_in_phase(lhs->get_genotype_call(sample), rhs->get_genotype_call(sample));
}

} // namespacedetail

template <typename Iterator>
void sort_genotype_alleles_by_phase_set(Iterator first_call, Iterator last_call, const SampleName& sample, const ReferenceGenome& reference)
{
    while (first_call != last_call) {
        const auto phased = [&] (const CallWrapper& call) { return detail::are_in_phase(call, *first_call, sample); };
        const auto last_phased_call = std::find_if_not(std::next(first_call), last_call, phased);
        sort_genotype_alleles(first_call, last_phased_call, sample, reference);
        first_call = last_phased_call;
    }
}

template <typename Range>
void sort_genotype_alleles(Range& calls, const SampleName& sample, const ReferenceGenome& reference)
{
    sort_genotype_alleles(std::begin(calls), std::end(calls), sample, reference);
}

template <typename Range>
void sort_genotype_alleles_by_phase_set(Range& calls, const SampleName& sample, const ReferenceGenome& reference)
{
    sort_genotype_alleles_by_phase_set(std::begin(calls), std::end(calls), sample, reference);
}

template <typename Iterator>
void sort_genotype_alleles(Iterator first_call, Iterator last_call, const std::vector<SampleName>& samples, const ReferenceGenome& reference)
{
    for (const auto& sample : samples) {
        sort_genotype_alleles(first_call, last_call, sample, reference);
    }
}

template <typename Iterator>
void sort_genotype_alleles_by_phase_set(Iterator first_call, Iterator last_call, const std::vector<SampleName>& samples, const ReferenceGenome& reference)
{
    for (const auto& sample : samples) {
        sort_genotype_alleles_by_phase_set(first_call, last_call, sample, reference);
    }
}

template <typename Range>
void sort_genotype_alleles(Range& calls, const std::vector<SampleName>& samples, const ReferenceGenome& reference)
{
    sort_genotype_alleles(std::begin(calls), std::end(calls), samples, reference);
}

template <typename Range>
void sort_genotype_alleles_by_phase_set(Range& calls, const std::vector<SampleName>& samples, const ReferenceGenome& reference)
{
    sort_genotype_alleles_by_phase_set(std::begin(calls), std::end(calls), samples, reference);
}

} // namespace octopus

#endif
