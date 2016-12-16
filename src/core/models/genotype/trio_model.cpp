// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "trio_model.hpp"

#include <iterator>
#include <algorithm>
#include <cmath>
#include <random>
#include <utility>
#include <cassert>
#include <string>
#include <iostream>

#include "utils/maths.hpp"
#include "germline_likelihood_model.hpp"

namespace octopus { namespace model {

TrioModel::TrioModel(const Trio& trio,
                     const PopulationPriorModel& prior_model,
                     const DeNovoModel& mutation_model,
                     Options options,
                     boost::optional<logging::DebugLogger> debug_log)
: trio_ {trio}
, prior_model_ {prior_model}
, mutation_model_ {mutation_model}
, options_ {options}
, debug_log_ {debug_log}
{}

namespace {

template <typename Container>
void clear(Container& c)
{
    c.clear();
    c.shrink_to_fit();
}

using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;

bool operator==(const GenotypeReference lhs, const GenotypeReference rhs)
{
    return lhs.get() == rhs.get();
}

bool operator<(const GenotypeReference lhs, const GenotypeReference rhs)
{
    return GenotypeLess{}(lhs.get(), rhs.get());
}

struct GenotypeRefProbabilityPair
{
    GenotypeReference genotype;
    double probability;
};

bool operator==(const GenotypeRefProbabilityPair& lhs, const GenotypeRefProbabilityPair& rhs) noexcept
{
    return lhs.probability == rhs.probability;
}
bool operator!=(const GenotypeRefProbabilityPair& lhs, const GenotypeRefProbabilityPair& rhs) noexcept
{
    return lhs.probability != rhs.probability;
}
bool operator<(const GenotypeRefProbabilityPair& lhs, const GenotypeRefProbabilityPair& rhs) noexcept
{
    return lhs.probability < rhs.probability;
}
bool operator>(const GenotypeRefProbabilityPair& lhs, const GenotypeRefProbabilityPair& rhs) noexcept
{
    return lhs.probability > rhs.probability;
}

struct ParentsProbabilityPair
{
    GenotypeReference maternal, paternal;
    double probability;
};

bool operator==(const ParentsProbabilityPair& lhs, const ParentsProbabilityPair& rhs) noexcept
{
    return lhs.probability == rhs.probability;
}
bool operator!=(const ParentsProbabilityPair& lhs, const ParentsProbabilityPair& rhs) noexcept
{
    return lhs.probability != rhs.probability;
}
bool operator<(const ParentsProbabilityPair& lhs, const ParentsProbabilityPair& rhs) noexcept
{
    return lhs.probability < rhs.probability;
}
bool operator>(const ParentsProbabilityPair& lhs, const ParentsProbabilityPair& rhs) noexcept
{
    return lhs.probability > rhs.probability;
}

auto compute_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                         const GermlineLikelihoodModel& model)
{
    std::vector<GenotypeRefProbabilityPair> result {};
    result.reserve(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::back_inserter(result),
                   [&model] (const auto& genotype) {
                       return GenotypeRefProbabilityPair {genotype, model.evaluate(genotype)};
                   });
    return result;
}

template <typename T>
bool all_equal(const std::vector<T>& v)
{
    return std::adjacent_find(std::cbegin(v), std::cend(v), std::not_equal_to<> {}) == std::cend(v);
}

template <typename Iterator>
double compute_removed_posterior_mass(Iterator first, Iterator first_removed, Iterator last)
{
    thread_local std::vector<double> likelihoods {};
    likelihoods.resize(std::distance(first, last));
    std::transform(first, last, std::begin(likelihoods),
                   [] (const auto& p) { return p.probability; });
    maths::normalise_exp(likelihoods);
    return std::accumulate(std::next(std::cbegin(likelihoods), std::distance(first, first_removed)),
                           std::cend(likelihoods), 0.0);
}

template <typename T>
auto compute_first_removable(const std::vector<T>& zipped, const double max_removed_mass)
{
    thread_local std::vector<double> likelihoods {};
    likelihoods.resize(zipped.size());
    std::transform(std::cbegin(zipped), std::cend(zipped), std::begin(likelihoods),
                   [] (const auto& p) { return p.probability; });
    maths::normalise_exp(likelihoods);
    // Do things in reverse order as floating points have better precision near zero
    std::partial_sum(std::crbegin(likelihoods), std::crend(likelihoods),
                     std::rbegin(likelihoods));
    const auto iter = std::upper_bound(std::crbegin(likelihoods), std::crend(likelihoods),
                                       max_removed_mass);
    const auto num_to_keep = std::distance(iter, std::crend(likelihoods));
    assert(num_to_keep <= zipped.size());
    return std::next(std::cbegin(zipped), num_to_keep);
}

template <typename T>
bool reduce(std::vector<T>& zipped, const std::size_t min_to_keep,
            const std::size_t max_to_keep, const double max_removed_mass)
{
    assert(min_to_keep <= max_to_keep);
    assert(0.0 <= max_removed_mass && max_removed_mass <= 1.0);
    bool hit_max_bound {false};
    if (zipped.size() <= min_to_keep) return hit_max_bound;
    if (all_equal(zipped)) {
        static std::default_random_engine gen {};
        std::shuffle(std::begin(zipped), std::end(zipped), gen);
        const auto min_can_keep = static_cast<std::size_t>(std::floor(zipped.size() * max_removed_mass));
        auto num_to_keep = std::max(min_to_keep, min_can_keep);
        if (num_to_keep > max_to_keep) {
            num_to_keep = max_to_keep;
            hit_max_bound = true;
        }
        assert(num_to_keep <= zipped.size());
        zipped.erase(std::next(std::cbegin(zipped), num_to_keep), std::cend(zipped));
    } else {
        const auto first_removed = std::next(std::begin(zipped), std::min(zipped.size(), min_to_keep));
        std::partial_sort(std::begin(zipped), first_removed, std::end(zipped), std::greater<> {});
        const auto removed_mass = compute_removed_posterior_mass(std::begin(zipped), first_removed, std::end(zipped));
        if (removed_mass > max_removed_mass) {
            std::sort(first_removed, std::end(zipped), std::greater<> {});
            zipped.erase(compute_first_removable(zipped, max_removed_mass), std::cend(zipped));
            if (zipped.size() > max_to_keep) {
                zipped.erase(std::next(std::cbegin(zipped), max_to_keep), std::cend(zipped));
                hit_max_bound = true;
            }
        } else {
            zipped.erase(first_removed, std::end(zipped));
        }
    }
    return hit_max_bound;
}

template <typename T>
auto reduce(std::vector<T>& zipped, const TrioModel::Options& options)
{
    return reduce(zipped, options.min_to_keep, options.max_to_keep,
                  options.max_removal_posterior_mass);
}

double probability_of_parents(const Genotype<Haplotype>& mother,
                              const Genotype<Haplotype>& father,
                              const PopulationPriorModel& model)
{
    thread_local std::vector<std::reference_wrapper<const Genotype<Haplotype>>> parental_genotypes {};
    if (parental_genotypes.empty()) {
        parental_genotypes.reserve(2);
        parental_genotypes.emplace_back(mother);
        parental_genotypes.emplace_back(father);
    } else {
        assert(parental_genotypes.size() == 2);
        parental_genotypes.front() = mother;
        parental_genotypes.back()  = father;
    }
    return model.evaluate(parental_genotypes);
}

auto join(const std::vector<GenotypeRefProbabilityPair>& maternal,
          const std::vector<GenotypeRefProbabilityPair>& paternal,
          const PopulationPriorModel& model)
{
    std::vector<ParentsProbabilityPair> result {};
    result.reserve(maternal.size() * paternal.size());
    for (const auto& m : maternal) {
        for (const auto& p : paternal) {
            result.push_back({m.genotype, p.genotype,
                                m.probability + p.probability
                                + probability_of_parents(m.genotype, p.genotype, model)});
        }
    }
    return result;
}

bool all_diploid(const Genotype<Haplotype>& child,
                 const Genotype<Haplotype>& mother,
                 const Genotype<Haplotype>& father)
{
    return is_diploid(child) && is_diploid(mother) && is_diploid(father);
}

double probability_of_child_given_parent(const Haplotype& child,
                                         const Genotype<Haplotype>& parent,
                                         const DeNovoModel& mutation_model)
{
    static const double ln2 {std::log(2)};
    const auto p1 = mutation_model.evaluate(child, parent[0]);
    const auto p2 = mutation_model.evaluate(child, parent[1]);
    return maths::log_sum_exp(p1, p2) - ln2;
}

double probability_of_child_given_parents(const Haplotype& child_from_mother,
                                          const Haplotype& child_from_father,
                                          const Genotype<Haplotype>& mother,
                                          const Genotype<Haplotype>& father,
                                          const DeNovoModel& mutation_model)
{
    return probability_of_child_given_parent(child_from_mother, mother, mutation_model)
            + probability_of_child_given_parent(child_from_father, father, mutation_model);
}

double probability_of_child_given_parents(const Genotype<Haplotype>& child,
                                          const Genotype<Haplotype>& mother,
                                          const Genotype<Haplotype>& father,
                                          const DeNovoModel& mutation_model)
{
    if (all_diploid(child, mother, father)) {
        static const double ln2 {std::log(2)};
        const auto p1 = probability_of_child_given_parents(child[0], child[1], mother, father, mutation_model);
        const auto p2 = probability_of_child_given_parents(child[1], child[0], mother, father, mutation_model);
        return maths::log_sum_exp(p1, p2) - ln2;
    }
    return 0; // TODO
}

using JointProbability = TrioModel::Latents::JointProbability;

auto joint_probability(const ParentsProbabilityPair& parents,
                       const GenotypeRefProbabilityPair& child,
                       const DeNovoModel& mutation_model)
{
    return parents.probability + child.probability
           + probability_of_child_given_parents(child.genotype, parents.maternal, parents.paternal, mutation_model);
}

auto join(const std::vector<ParentsProbabilityPair>& parents,
          const std::vector<GenotypeRefProbabilityPair>& child,
          const DeNovoModel& mutation_model)
{
    std::vector<JointProbability> result {};
    result.reserve(parents.size() * child.size());
    for (const auto& p : parents) {
        for (const auto& c : child) {
            result.push_back({p.maternal, p.paternal, c.genotype, joint_probability(p, c, mutation_model)});
        }
    }
    return result;
}

auto extract_probabilities(const std::vector<JointProbability>& joint_likelihoods)
{
    std::vector<double> result(joint_likelihoods.size());
    std::transform(std::cbegin(joint_likelihoods), std::cend(joint_likelihoods),
                   std::begin(result), [] (const auto& p) { return p.probability; });
    return result;
}

auto normalise_exp(std::vector<JointProbability>& joint_likelihoods)
{
    auto likelihoods = extract_probabilities(joint_likelihoods);
    const auto norm = maths::normalise_exp(likelihoods);
    auto iter = std::cbegin(likelihoods);
    for (auto& p : joint_likelihoods) p.probability = *iter++;
    return norm;
}

} // namespace

namespace debug {

template <typename S>
void print(S&& stream, const std::string& sample, std::vector<GenotypeRefProbabilityPair> ps,
           std::size_t n = 5);
void print(const std::string& sample, const std::vector<GenotypeRefProbabilityPair>& ps,
           std::size_t n = 5);
template <typename S>
void print(S&& stream, std::vector<ParentsProbabilityPair> ps, std::size_t n = 5);
void print(std::vector<ParentsProbabilityPair> ps, std::size_t n = 5);
template <typename S>
void print(S&& stream, std::vector<JointProbability> ps, std::size_t n = 5);
void print(std::vector<JointProbability> ps, std::size_t n = 5);

} // namespace debug

TrioModel::InferredLatents
TrioModel::evaluate(const GenotypeVector& maternal_genotypes,
                    const GenotypeVector& paternal_genotypes,
                    const GenotypeVector& child_genotypes,
                    const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!maternal_genotypes.empty() && !paternal_genotypes.empty() && !child_genotypes.empty());
    const GermlineLikelihoodModel likelihood_model {haplotype_likelihoods};
    bool overflowed {false};
    haplotype_likelihoods.prime(trio_.mother());
    auto maternal_likelihoods = compute_likelihoods(maternal_genotypes, likelihood_model);
    overflowed |= reduce(maternal_likelihoods, options_);
    assert(!maternal_likelihoods.empty());
    haplotype_likelihoods.prime(trio_.father());
    auto paternal_likelihoods = compute_likelihoods(paternal_genotypes, likelihood_model);
    overflowed |= reduce(paternal_likelihoods, options_);
    assert(!paternal_likelihoods.empty());
    if (debug_log_) {
        debug::print(stream(*debug_log_), "maternal", maternal_likelihoods);
        debug::print(stream(*debug_log_), "paternal", paternal_likelihoods);
    }
    auto parents_joint_likelihoods = join(maternal_likelihoods, paternal_likelihoods, prior_model_);
    overflowed |= reduce(parents_joint_likelihoods, options_);
    assert(!parents_joint_likelihoods.empty());
    clear(maternal_likelihoods);
    clear(paternal_likelihoods);
    if (debug_log_) debug::print(stream(*debug_log_), parents_joint_likelihoods);
    haplotype_likelihoods.prime(trio_.child());
    auto child_likelihoods = compute_likelihoods(child_genotypes, likelihood_model);
    overflowed |= reduce(child_likelihoods, options_);
    assert(!child_likelihoods.empty());
    if (debug_log_) debug::print(stream(*debug_log_), "child", child_likelihoods);
    auto joint_likelihoods = join(parents_joint_likelihoods, child_likelihoods, mutation_model_);
    clear(parents_joint_likelihoods);
    clear(child_likelihoods);
    const auto evidence = normalise_exp(joint_likelihoods);
    if (debug_log_) debug::print(stream(*debug_log_), joint_likelihoods);
    return {std::move(joint_likelihoods), evidence, overflowed};
}

namespace debug {

template <typename S>
void print(S&& stream, const std::string& sample, std::vector<GenotypeRefProbabilityPair> ps,
           const std::size_t n)
{
    stream << sample << ":\n";
    const auto nth = std::next(std::begin(ps), std::min(ps.size(), n));
    std::partial_sort(std::begin(ps), nth, std::end(ps),
                      [] (const auto& lhs, const auto& rhs) {
                          return lhs.probability > rhs.probability;
                      });
    std::for_each(std::begin(ps), nth, [&] (const auto& p) {
        using octopus::debug::print_variant_alleles;
        print_variant_alleles(stream, p.genotype);
        stream << " " << p.probability << "\n";
    });
}

void print(const std::string& sample, const std::vector<GenotypeRefProbabilityPair>& ps,
           std::size_t n)
{
    print(std::cout, sample, std::move(ps), n);
}

template <typename S>
void print(S&& stream, std::vector<ParentsProbabilityPair> ps, const std::size_t n)
{
    stream << "joint parents" << ":\n";
    const auto nth = std::next(std::begin(ps), std::min(ps.size(), n));
    std::partial_sort(std::begin(ps), nth, std::end(ps),
                      [] (const auto& lhs, const auto& rhs) {
                          return lhs.probability > rhs.probability;
                      });
    std::for_each(std::begin(ps), nth, [&] (const auto& p) {
        using octopus::debug::print_variant_alleles;
        print_variant_alleles(stream, p.maternal);
        stream << " | ";
        print_variant_alleles(stream, p.paternal);
        stream << " " << p.probability << "\n";
    });
}

void print(std::vector<ParentsProbabilityPair> ps, std::size_t n)
{
    print(std::cout, std::move(ps), n);
}

template <typename S>
void print(S&& stream, std::vector<JointProbability> ps, const std::size_t n)
{
    stream << "trio top (maternal | paternal | child):" << '\n';
    const auto nth = std::next(std::begin(ps), std::min(ps.size(), n));
    std::partial_sort(std::begin(ps), nth, std::end(ps),
                      [] (const auto& lhs, const auto& rhs) {
                          return lhs.probability > rhs.probability;
                      });
    std::for_each(std::begin(ps), nth, [&] (const auto& p) {
        using octopus::debug::print_variant_alleles;
        print_variant_alleles(stream, p.maternal);
        stream << " | ";
        print_variant_alleles(stream, p.paternal);
        stream << " | ";
        print_variant_alleles(stream, p.child);
        stream << " " << p.probability << "\n";
    });
}

void print(std::vector<JointProbability> ps, std::size_t n)
{
    print(std::cout, std::move(ps), n);
}

} // namespace debug

} // namespace model
} // namespace octopus
