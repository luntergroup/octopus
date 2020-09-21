// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "hardy_weinberg_model.hpp"

#include <utility>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <cmath>

#include "utils/maths.hpp"

#include <iostream>

namespace octopus {

HardyWeinbergModel::HardyWeinbergModel(IndexedHaplotype<> reference)
: reference_ {std::move(reference)}
, haplotype_frequencies_ {}
, empirical_ {false}
{}

HardyWeinbergModel::HardyWeinbergModel(HaplotypeFrequencyVector haplotype_frequencies)
: reference_ {}
, haplotype_frequencies_ {std::move(haplotype_frequencies)}
, empirical_ {true}
{}

void HardyWeinbergModel::set_frequencies(HaplotypeFrequencyVector haplotype_frequencies)
{
    haplotype_frequencies_ = std::move(haplotype_frequencies);
    empirical_ = true;
}

HardyWeinbergModel::HaplotypeFrequencyVector& HardyWeinbergModel::frequencies() noexcept
{
    return haplotype_frequencies_;
}

const HardyWeinbergModel::HaplotypeFrequencyVector& HardyWeinbergModel::frequencies() const noexcept
{
    return haplotype_frequencies_;
}

namespace {

auto ln_hardy_weinberg_haploid(const Genotype<IndexedHaplotype<>>& genotype,
                               const HardyWeinbergModel::HaplotypeFrequencyVector& haplotype_frequencies)
{
    return std::log(haplotype_frequencies[index_of(genotype[0])]);
}

auto ln_hardy_weinberg_diploid(const Genotype<IndexedHaplotype<>>& genotype,
                               const HardyWeinbergModel::HaplotypeFrequencyVector& haplotype_frequencies)
{
    if (genotype[0] == genotype[1]) {
        return 2 * std::log(haplotype_frequencies[index_of(genotype[0])]);
    }
    static const double ln2 {std::log(2.0)};
    return std::log(haplotype_frequencies[index_of(genotype[0])]) + std::log(haplotype_frequencies[index_of(genotype[1])]) + ln2;
}

auto ln_hardy_weinberg_polyploid(const Genotype<IndexedHaplotype<>>& genotype,
                                 const HardyWeinbergModel::HaplotypeFrequencyVector& haplotype_frequencies)
{
    std::vector<unsigned> counts(haplotype_frequencies.size());
    for (const auto& haplotype : genotype) ++counts[index_of(haplotype)];
    return maths::log_multinomial_pdf<>(counts, haplotype_frequencies);
}

template <typename Range>
auto sum(const Range& values) noexcept
{
    using T = typename Range::value_type;
    return std::accumulate(std::cbegin(values), std::cend(values), T {0});
}

template <typename T>
std::vector<T> to_frequencies(const std::vector<unsigned>& counts)
{
    std::vector<T> result(counts.size());
    const auto norm = static_cast<T>(sum(counts));
    std::transform(std::cbegin(counts), std::cend(counts), std::begin(result),
                   [norm] (auto count) noexcept { return static_cast<T>(count) / norm; });
    return result;
}

} // namespace

HardyWeinbergModel::LogProbability HardyWeinbergModel::evaluate(const Genotype<IndexedHaplotype<>>& genotype) const
{
    if (empirical_) {
        switch (genotype.ploidy()) {
            case 1 : return ln_hardy_weinberg_haploid(genotype, haplotype_frequencies_);
            case 2 : return ln_hardy_weinberg_diploid(genotype, haplotype_frequencies_);
            default: return ln_hardy_weinberg_polyploid(genotype, haplotype_frequencies_);
        }
    } else {
        static const LogProbability ln2 {std::log(2.0)}, ln3 {std::log(3.0)};
        if (is_haploid(genotype)) {
            return reference_ && contains(genotype, *reference_) ? -ln2 : 0.0;
        }
        if (is_diploid(genotype)) {
            if (reference_ && contains(genotype, *reference_)) {
                return is_homozygous(genotype) ? -ln2 : -ln3;
            } else {
                return is_homozygous(genotype) ? -ln2 : 0.0;
            }
        }
        auto counts = unique_counts(genotype);
        if (reference_ && !contains(genotype, *reference_)) {
            counts.push_back(1);
        }
        auto probs = to_frequencies<LogProbability>(counts);
        return maths::log_multinomial_pdf(counts, probs);
    }
}

namespace {

template <typename T> const T& get(const T& value) noexcept { return value; }
template <typename T> const T& get(std::reference_wrapper<const T> value) noexcept { return value.get(); }

template <typename Range>
auto sum_plodies(const Range& genotypes) noexcept
{
    return std::accumulate(std::cbegin(genotypes), std::cend(genotypes), 0u,
                           [] (auto curr, const auto& g) { return curr + get(g).ploidy(); });
}

template <typename Range>
void fill_frequencies(const Range& genotypes, HardyWeinbergModel::HaplotypeFrequencyVector& result)
{
    decltype(index_of(get(genotypes[0])[0])) max_haplotype_idx {0}, n {0};
    for (const auto& genotype : genotypes) {
        for (const auto& haplotype : get(genotype)) {
            max_haplotype_idx = std::max(max_haplotype_idx, index_of(haplotype));
            ++n;
        }
    }
    result.resize(max_haplotype_idx + 1);
    const auto weight = 1.0 / n;
    for (const auto& genotype : genotypes) {
        for (const auto& haplotype : get(genotype)) {
            result[index_of(haplotype)] += weight;
        }
    }
}

template <typename Range>
auto joint_evaluate(const Range& genotypes, const HardyWeinbergModel& model)
{
    return std::accumulate(std::cbegin(genotypes), std::cend(genotypes), 0.0,
                           [&] (auto curr, const auto& genotype) { return curr + model.evaluate(genotype); });
}

} // namespace

HardyWeinbergModel::LogProbability HardyWeinbergModel::evaluate(const std::vector<Genotype<IndexedHaplotype<>>>& genotypes) const
{
    if (empirical_) {
        return joint_evaluate(genotypes, *this);
    } else {
        fill_frequencies(genotypes, haplotype_frequencies_);
        empirical_ = true;
        auto result = joint_evaluate(genotypes, *this);
        haplotype_frequencies_.clear();
        empirical_ = false;
        return result;
    }
}

HardyWeinbergModel::LogProbability HardyWeinbergModel::evaluate(const std::vector<GenotypeReference>& genotypes) const
{
    if (empirical_) {
        return joint_evaluate(genotypes, *this);
    } else {
        fill_frequencies(genotypes, haplotype_frequencies_);
        empirical_ = true;
        auto result = joint_evaluate(genotypes, *this);
        haplotype_frequencies_.clear();
        empirical_ = false;
        return result;
    }
}

} // namespace
