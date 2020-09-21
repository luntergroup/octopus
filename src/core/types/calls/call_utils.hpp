// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef call_utils_hpp
#define call_utils_hpp

#include <vector>
#include <iterator>
#include <cstddef>
#include <algorithm>
#include <numeric>

#include "core/types/genotype.hpp"
#include "core/types/haplotype.hpp"
#include "io/reference/reference_genome.hpp"
#include "config/common.hpp"
#include "call.hpp"
#include "call_wrapper.hpp"

namespace octopus {

namespace detail {

template <typename Iterator>
std::size_t count_calls(const std::vector<std::pair<Iterator, Iterator>>& calls)
{
    const static auto add_count = [] (auto total, const auto& p) { return total + std::distance(p.first, p.second); };
    return std::accumulate(std::cbegin(calls), std::cend(calls), std::size_t {0}, add_count);
}

template <typename Iterator>
auto phase_set_ploidy(std::vector<std::pair<Iterator, Iterator>>& calls, const SampleName& sample)
{
    assert(!calls.empty());
    const CallWrapper& call {*calls.front().first};
    return call->get_genotype_call(sample).genotype.ploidy();
}

struct AlleleExtractor
{
    explicit AlleleExtractor(unsigned index) noexcept : index_ {index} {}
    const Allele& operator()(const Genotype<Allele>& genotype) const noexcept { return genotype[index_]; }
private:
    unsigned index_;
};

std::vector<unsigned>
compute_haplotype_order(std::vector<Genotype<Allele>>& genotypes, const ReferenceGenome& reference);

} // namespace detail

template <typename Iterator>
void sort_genotype_alleles(std::vector<std::pair<Iterator, Iterator>>& calls, const SampleName& sample, const ReferenceGenome& reference)
{
    const auto num_calls = detail::count_calls(calls);
    if (num_calls > 0 && detail::phase_set_ploidy(calls, sample) > 0) {
        std::vector<Genotype<Allele>> genotypes {};
        genotypes.reserve(num_calls);
        for (const auto& p : calls) {
            std::for_each(p.first, p.second, [&] (CallWrapper& call) {
                genotypes.push_back(call->get_genotype_call(sample).genotype);
            });
        }
        const auto new_allele_order = detail::compute_haplotype_order(genotypes, reference);
        if (!std::is_sorted(std::cbegin(new_allele_order), std::cend(new_allele_order))) {
            for (auto& p : calls) {
                std::for_each(p.first, p.second, [&] (CallWrapper& call) {
                    call->reorder_genotype(sample, new_allele_order);
                });
            }
        }
    }
}

namespace detail {

bool are_same_phase_set(const Call::GenotypeCall& lhs, const Call::GenotypeCall& rhs);
bool are_same_phase_set(const CallWrapper& lhs, const CallWrapper& rhs, const SampleName& sample);
bool have_different_phase_sets(const Call::GenotypeCall& lhs, const Call::GenotypeCall& rhs);
bool have_different_phase_sets(const CallWrapper& lhs, const CallWrapper& rhs, const SampleName& sample);

} // namespace detail

template <typename Iterator>
std::vector<std::vector<std::pair<Iterator, Iterator>>>
segment_by_phase_set(Iterator first_call, Iterator last_call, const SampleName& sample)
{
    std::vector<std::vector<std::pair<Iterator, Iterator>>> result {};
    std::vector<Iterator> phase_set_indicators {};
    while (first_call != last_call) {
        const auto is_phased = [&] (const CallWrapper& call) { return call->is_phased(sample); };
        const auto first_phased_call = std::find_if(first_call, last_call, is_phased);
        if (first_phased_call != last_call) {
            const auto different_phase_sets = [&] (const CallWrapper& call) { 
                return detail::have_different_phase_sets(*first_phased_call, call, sample); };
            const auto last_phased_call = std::find_if(std::next(first_phased_call), last_call, different_phase_sets);
            const auto same_phase_set = [&] (const auto& indicator_itr) {
                 return detail::are_same_phase_set(*indicator_itr, *first_phased_call, sample); };
            const auto indicator_itr = std::find_if(std::cbegin(phase_set_indicators), std::cend(phase_set_indicators), same_phase_set);
            if (indicator_itr == std::cend(phase_set_indicators)) {
                result.push_back({{first_call, last_phased_call}});
                phase_set_indicators.push_back(first_phased_call);
            } else {
                const auto result_idx = static_cast<std::size_t>(std::distance(std::cbegin(phase_set_indicators), indicator_itr));
                result[result_idx].push_back({first_call, last_phased_call});
            }
            first_call = last_phased_call;
        } else {
            result.push_back({{first_call, last_call}});
            break;
        }
    }
    return result;
}

template <typename Iterator>
void sort_genotype_alleles_by_phase_set(Iterator first_call, Iterator last_call, const SampleName& sample, const ReferenceGenome& reference)
{
    auto phase_sets = segment_by_phase_set(first_call, last_call, sample);
    for (auto& phase_set : phase_sets) {
        sort_genotype_alleles(phase_set, sample, reference);
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
