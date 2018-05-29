// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "hardy_weinberg_model.hpp"

#import <utility>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <cmath>

#include "utils/maths.hpp"

#include <iostream>

namespace octopus {

HardyWeinbergModel::HardyWeinbergModel(Haplotype reference)
: reference_ {std::move(reference)}
, reference_idx_ {}
, haplotype_frequencies_ {}
, haplotype_idx_frequencies_ {}
, empirical_ {false}
{}

HardyWeinbergModel::HardyWeinbergModel(unsigned reference_idx)
: reference_ {}
, reference_idx_ {reference_idx}
, haplotype_frequencies_ {}
, haplotype_idx_frequencies_ {}
, empirical_ {false}
{}

HardyWeinbergModel::HardyWeinbergModel(HaplotypeFrequencyMap haplotype_frequencies)
: reference_ {}
, reference_idx_ {}
, haplotype_frequencies_ {std::move(haplotype_frequencies)}
, haplotype_idx_frequencies_ {}
, empirical_ {true}
{}

HardyWeinbergModel::HardyWeinbergModel(HaplotypeFrequencyVector haplotype_frequencies)
: reference_ {}
, reference_idx_ {}
, haplotype_frequencies_ {}
, haplotype_idx_frequencies_ {std::move(haplotype_frequencies)}
, empirical_ {true}
{}

void HardyWeinbergModel::set_frequencies(HaplotypeFrequencyMap haplotype_frequencies)
{
    haplotype_frequencies_ = std::move(haplotype_frequencies);
    empirical_ = true;
}

void HardyWeinbergModel::set_frequencies(HaplotypeFrequencyVector haplotype_frequencies)
{
    haplotype_idx_frequencies_ = std::move(haplotype_frequencies);
    empirical_ = true;
}

HardyWeinbergModel::HaplotypeFrequencyMap& HardyWeinbergModel::frequencies() noexcept
{
    return haplotype_frequencies_;
}

HardyWeinbergModel::HaplotypeFrequencyVector& HardyWeinbergModel::index_frequencies() noexcept
{
    return haplotype_idx_frequencies_;
}

namespace {

auto ln_hardy_weinberg_haploid(const Genotype<Haplotype>& genotype,
                               const HardyWeinbergModel::HaplotypeFrequencyMap& haplotype_frequencies)
{
    return std::log(haplotype_frequencies.at(genotype[0]));
}

auto ln_hardy_weinberg_diploid(const Genotype<Haplotype>& genotype,
                               const HardyWeinbergModel::HaplotypeFrequencyMap& haplotype_frequencies)
{
    if (genotype.is_homozygous()) {
        return 2 * std::log(haplotype_frequencies.at(genotype[0]));
    }
    static const double ln2 {std::log(2.0)};
    return std::log(haplotype_frequencies.at(genotype[0])) + std::log(haplotype_frequencies.at(genotype[1])) + ln2;
}

auto ln_hardy_weinberg_polyploid(const Genotype<Haplotype>& genotype,
                                 const HardyWeinbergModel::HaplotypeFrequencyMap& haplotype_frequencies)
{
    auto unique_haplotypes = genotype.copy_unique();
    std::vector<unsigned> occurences {};
    occurences.reserve(unique_haplotypes.size());
    double r {0};
    for (const auto& haplotype : unique_haplotypes) {
        auto num_occurences = genotype.count(haplotype);
        occurences.push_back(num_occurences);
        r += num_occurences * std::log(haplotype_frequencies.at(haplotype));
    }
    return maths::log_multinomial_coefficient<double>(occurences) + r;
}

auto ln_hardy_weinberg_haploid(const GenotypeIndex& genotype,
                               const HardyWeinbergModel::HaplotypeFrequencyVector& haplotype_frequencies)
{
    return std::log(haplotype_frequencies[genotype[0]]);
}

auto ln_hardy_weinberg_diploid(const GenotypeIndex& genotype,
                               const HardyWeinbergModel::HaplotypeFrequencyVector& haplotype_frequencies)
{
    if (genotype[0] == genotype[1]) {
        return 2 * std::log(haplotype_frequencies[genotype[0]]);
    }
    static const double ln2 {std::log(2.0)};
    return std::log(haplotype_frequencies[genotype[0]]) + std::log(haplotype_frequencies[genotype[1]]) + ln2;
}

auto ln_hardy_weinberg_polyploid(const GenotypeIndex& genotype,
                                 const HardyWeinbergModel::HaplotypeFrequencyVector& haplotype_frequencies)
{
    std::vector<unsigned> counts(haplotype_frequencies.size());
    for (auto idx : genotype) ++counts[idx];
    return maths::log_multinomial_pdf<>(counts, haplotype_frequencies);
}

template <typename Range>
auto sum(const Range& values) noexcept
{
    using T = typename Range::value_type;
    return std::accumulate(std::cbegin(values), std::cend(values), T {0});
}

std::vector<double> to_frequencies(const std::vector<unsigned>& counts)
{
    std::vector<double> result(counts.size());
    const auto norm = static_cast<double>(sum(counts));
    std::transform(std::cbegin(counts), std::cend(counts), std::begin(result),
                   [norm] (auto count) noexcept { return static_cast<double>(count) / norm; });
    return result;
}

} // namespace

double HardyWeinbergModel::evaluate(const Genotype<Haplotype>& genotype) const
{
    if (empirical_) {
        switch (genotype.ploidy()) {
            case 1 : return ln_hardy_weinberg_haploid(genotype, haplotype_frequencies_);
            case 2 : return ln_hardy_weinberg_diploid(genotype, haplotype_frequencies_);
            default: return ln_hardy_weinberg_polyploid(genotype, haplotype_frequencies_);
        }
    } else {
        static const double ln2 {std::log(2.0)}, ln3 {std::log(3.0)};
        if (is_haploid(genotype)) {
            return reference_ && genotype.contains(*reference_) ? -ln2 : 0.0;
        }
        if (is_diploid(genotype)) {
            if (reference_ && genotype.contains(*reference_)) {
                return genotype.is_homozygous() ? -ln2 : -ln3;
            } else {
                return genotype.is_homozygous() ? -ln2 : 0.0;
            }
        }
        auto counts = genotype.unique_counts();
        if (reference_ && !genotype.contains(*reference_)) {
            counts.push_back(1);
        }
        auto probs = to_frequencies(counts);
        return maths::log_multinomial_pdf(counts, probs);
    }
}

namespace {

template <typename Range>
void unique_counts(const Range& range, std::vector<unsigned>& result)
{
    for (auto itr = std::cbegin(range), last = std::cend(range); itr != last;) {
        auto next = std::find_if_not(std::next(itr), last, [itr] (const auto& x) { return x == *itr; });
        result.push_back(std::distance(itr, next));
        itr = next;
    }
}

} // namespace

double HardyWeinbergModel::evaluate(const GenotypeIndex& genotype) const
{
    assert(!genotype.empty());
    if (empirical_) {
        switch (genotype.size()) {
            case 1 : return ln_hardy_weinberg_haploid(genotype, haplotype_idx_frequencies_);
            case 2 : return ln_hardy_weinberg_diploid(genotype, haplotype_idx_frequencies_);
            default: return ln_hardy_weinberg_polyploid(genotype, haplotype_idx_frequencies_);
        }
    } else {
        static const double ln2 {std::log(2.0)}, ln3 {std::log(3.0)};
        if (genotype.size() == 1) {
            return reference_idx_ && genotype[0] == *reference_idx_ ? -ln2 : 0.0;
        }
        if (genotype.size() == 2) {
            if (reference_idx_ && !(genotype[0] == *reference_idx_ || genotype[1] == *reference_idx_)) {
                return genotype[0] == genotype[1] ? -ln2 : -ln3;
            } else {
                return genotype[0] == genotype[1] ? -ln2 : 0.0;
            }
        }
        std::vector<unsigned> counts {};
        counts.reserve(genotype.size());
        if (std::is_sorted(std::cbegin(genotype), std::cend(genotype))) {
            unique_counts(genotype, counts);
        } else {
            auto sorted_genotype = genotype;
            std::sort(std::begin(sorted_genotype), std::end(sorted_genotype));
            unique_counts(sorted_genotype, counts);
        }
        if (reference_idx_ && std::find(std::cbegin(genotype), std::cend(genotype), *reference_idx_) == std::cend(genotype)) {
            counts.push_back(1);
        }
        auto probs = to_frequencies(counts);
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
void fill_frequencies(const Range& genotypes, HardyWeinbergModel::HaplotypeFrequencyMap& result)
{
    const auto n = sum_plodies(genotypes);
    const auto weight = 1.0 / n;
    for (const auto& genotype : genotypes) {
        for (const auto& haplotype : get(genotype)) {
            result[haplotype] += weight;
        }
    }
}

template <typename Range>
void fill_frequencies(const Range& genotypes, HardyWeinbergModel::HaplotypeFrequencyVector& result)
{
    unsigned max_haplotype_idx {0}, n {0};
    for (const auto& genotype : genotypes) {
        for (auto haplotype_idx : get(genotype)) {
            max_haplotype_idx = std::max(max_haplotype_idx, haplotype_idx);
            ++n;
        }
    }
    result.resize(max_haplotype_idx + 1);
    const auto weight = 1.0 / n;
    for (const auto& genotype : genotypes) {
        for (auto haplotype_idx : get(genotype)) {
            result[haplotype_idx] += weight;
        }
    }
}

template <typename Range>
double joint_evaluate(const Range& genotypes, const HardyWeinbergModel& model)
{
    return std::accumulate(std::cbegin(genotypes), std::cend(genotypes), 0.0,
                           [&model] (auto curr, const auto& genotype) { return curr + model.evaluate(get(genotype)); });
}

} // namespace

double HardyWeinbergModel::evaluate(const std::vector<Genotype<Haplotype>>& genotypes) const
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

double HardyWeinbergModel::evaluate(const GenotypeReferenceVector& genotypes) const
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

double HardyWeinbergModel::evaluate(const GenotypeIndexVector& genotypes) const
{
    if (empirical_) {
        return joint_evaluate(genotypes, *this);
    } else {
        fill_frequencies(genotypes, haplotype_idx_frequencies_);
        empirical_ = true;
        auto result = joint_evaluate(genotypes, *this);
        haplotype_idx_frequencies_.clear();
        empirical_ = false;
        return result;
    }
}

double HardyWeinbergModel::evaluate(const GenotypeIndexReferenceVector& genotypes) const
{
    if (empirical_) {
        return joint_evaluate(genotypes, *this);
    } else {
        fill_frequencies(genotypes, haplotype_idx_frequencies_);
        empirical_ = true;
        auto result = joint_evaluate(genotypes, *this);
        haplotype_idx_frequencies_.clear();
        empirical_ = false;
        return result;
    }
}

} // namespace
