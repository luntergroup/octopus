// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef call_utils_hpp
#define call_utils_hpp

#include <vector>
#include <iterator>
#include <cstddef>
#include <algorithm>

#include "core/types/genotype.hpp"
#include "config/common.hpp"
#include "call_wrapper.hpp"

namespace octopus {

template <typename Iterator>
void sort_genotype_alleles(Iterator first_call, Iterator last_call, const SampleName& sample)
{
    const auto num_calls = static_cast<std::size_t>(std::distance(first_call, last_call));
    if (num_calls > 0) {
        std::vector<Genotype<Allele>> genotypes {};
        genotypes.reserve(num_calls);
        std::for_each(first_call, last_call, [&] (CallWrapper& call) {
            genotypes.push_back(call->get_genotype_call(sample).genotype);
        });
        sort_alleles_in_haplotype_order(genotypes);
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
void sort_genotype_alleles_by_phase_set(Iterator first_call, Iterator last_call, const SampleName& sample)
{
    while (first_call != last_call) {
        const auto phased = [&] (const CallWrapper& call) { return detail::are_in_phase(call, *first_call, sample); };
        const auto last_phased_call = std::find_if_not(std::next(first_call), last_call, phased);
        sort_genotype_alleles(first_call, last_phased_call, sample);
        first_call = last_phased_call;
    }
}

template <typename Range>
void sort_genotype_alleles(Range& calls, const SampleName& sample)
{
    sort_genotype_alleles(std::begin(calls), std::end(calls), sample);
}

template <typename Range>
void sort_genotype_alleles_by_phase_set(Range& calls, const SampleName& sample)
{
    sort_genotype_alleles_by_phase_set(std::begin(calls), std::end(calls), sample);
}

template <typename Iterator>
void sort_genotype_alleles(Iterator first_call, Iterator last_call, const std::vector<SampleName>& samples)
{
    for (const auto& sample : samples) {
        sort_genotype_alleles(first_call, last_call, sample);
    }
}

template <typename Iterator>
void sort_genotype_alleles_by_phase_set(Iterator first_call, Iterator last_call, const std::vector<SampleName>& samples)
{
    for (const auto& sample : samples) {
        sort_genotype_alleles_by_phase_set(first_call, last_call, sample);
    }
}

template <typename Range>
void sort_genotype_alleles(Range& calls, const std::vector<SampleName>& samples)
{
    sort_genotype_alleles(std::begin(calls), std::end(calls), samples);
}

template <typename Range>
void sort_genotype_alleles_by_phase_set(Range& calls, const std::vector<SampleName>& samples)
{
    sort_genotype_alleles_by_phase_set(std::begin(calls), std::end(calls), samples);
}

} // namespace octopus

#endif
