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

struct GenotypeRefProbabilityPairGenotypeEqual
{
    bool operator()(const GenotypeRefProbabilityPair& lhs, const GenotypeRefProbabilityPair& rhs) const
    {
        return lhs.genotype == rhs.genotype;
    }
};

struct GenotypeRefProbabilityPairGenotypeLess
{
    bool operator()(const GenotypeRefProbabilityPair& lhs, const GenotypeRefProbabilityPair& rhs) const
    {
        return lhs.genotype < rhs.genotype;
    }
};

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

auto compute_posteriors(const std::vector<GenotypeRefProbabilityPair>& likelihoods,
                        const PopulationPriorModel& prior_model)
{
    std::vector<double> posteriors(likelihoods.size());
    std::transform(std::cbegin(likelihoods), std::cend(likelihoods), std::begin(posteriors),
                   [&] (const auto& p) {
                       return p.probability + prior_model.evaluate(std::vector<GenotypeReference> {p.genotype});
                   });
    maths::normalise_exp(posteriors);
    std::vector<GenotypeRefProbabilityPair> result {};
    result.reserve(posteriors.size());
    std::transform(std::cbegin(likelihoods), std::cend(likelihoods), std::cbegin(posteriors), std::back_inserter(result),
                   [] (const auto& p, const auto posterior) noexcept {
                       return GenotypeRefProbabilityPair {p.genotype, posterior};
                   });
    return result;
}

template <typename T>
bool all_equal(const std::vector<T>& v)
{
    return std::adjacent_find(std::cbegin(v), std::cend(v), std::not_equal_to<> {}) == std::cend(v);
}

template <typename Iterator>
double compute_lost_posterior_mass(Iterator first, Iterator first_removed, Iterator last)
{
    thread_local std::vector<double> likelihoods {};
    likelihoods.resize(std::distance(first, last));
    std::transform(first, last, std::begin(likelihoods), [] (const auto& p) { return p.probability; });
    maths::normalise_exp(likelihoods);
    return std::accumulate(std::cbegin(likelihoods), std::next(std::cbegin(likelihoods), std::distance(first, first_removed)), 0.0);
}

template <typename T>
auto compute_min_joint(const std::vector<T>& zipped, const double max_mass_loss)
{
    assert(!zipped.empty());
    std::vector<double> posteriors {};
    posteriors.resize(zipped.size());
    std::transform(std::cbegin(zipped), std::cend(zipped), std::begin(posteriors),
                   [] (const auto& p) { return p.probability; });
    maths::normalise_exp(posteriors);
    // Do things in reverse order as floating points have better precision near zero
    std::partial_sum(std::crbegin(posteriors), std::crend(posteriors), std::rbegin(posteriors));
    const auto iter = std::upper_bound(std::crbegin(posteriors), std::crend(posteriors), max_mass_loss);
    const auto num_to_keep = static_cast<std::size_t>(std::distance(iter, std::crend(posteriors)));
    assert(num_to_keep <= zipped.size());
    return std::max(num_to_keep, std::size_t {1});
}

template <typename T>
typename std::vector<T>::iterator
reduce(std::vector<T>& zipped, std::size_t max_joint, const double max_mass_loss)
{
    max_joint = std::min(max_joint, zipped.size());
    if (all_equal(zipped)) {
        static std::default_random_engine gen {};
        std::shuffle(std::begin(zipped), std::end(zipped), gen);
        const auto min_joint = static_cast<std::size_t>(std::floor(zipped.size() * (1 - max_mass_loss)));
        return std::next(std::begin(zipped), std::min(max_joint, min_joint));
    } else {
        auto result = std::next(std::begin(zipped), max_joint);
        std::partial_sort(std::begin(zipped), result, std::end(zipped), std::greater<> {});
        const auto lost_mass = compute_lost_posterior_mass(std::begin(zipped), result, std::end(zipped));
        if (lost_mass > max_mass_loss) {
            const auto min_joint = compute_min_joint(zipped, max_mass_loss);
            if (min_joint < max_joint) {
                result = std::next(std::begin(zipped), min_joint);
            }
        }
        assert(result <= std::end(zipped));
        return result;
    }
}

unsigned get_sample_reduction_count(const std::size_t n)
{
    return static_cast<unsigned>(std::sqrt(n));
}

template <typename T>
struct ReducedVectorMap
{
    using Iterator = typename std::vector<T>::const_iterator;
    Iterator first, last_to_partially_join, last_to_join, last;
};

template <typename T>
auto make_reduction_map(const std::vector<T>& elements,
                        const typename std::vector<T>::const_iterator last_to_combine,
                        const TrioModel::Options& options)
{
    if (std::distance(std::cbegin(elements), last_to_combine) > 1) {
        return ReducedVectorMap<T> {std::cbegin(elements), std::next(std::cbegin(elements), 1), last_to_combine, std::cend(elements)};
    } else {
        return ReducedVectorMap<T> {std::cbegin(elements), last_to_combine, last_to_combine, std::cend(elements)};
    }
}

template <typename T>
auto reduce(std::vector<T>& zipped, const TrioModel::Options& options)
{
    auto last_to_join = reduce(zipped,
                               get_sample_reduction_count(options.max_joint_genotypes),
                               options.max_mass_loss);
    return make_reduction_map(zipped, last_to_join, options);
}

auto reduce(std::vector<ParentsProbabilityPair>& zipped, const TrioModel::Options& options)
{
    auto last_to_join = reduce(zipped,
                               get_sample_reduction_count(options.max_joint_genotypes),
                               std::pow(options.max_mass_loss, 2));
    return make_reduction_map(zipped, last_to_join, options);
}

bool are_same_ploidy(const std::vector<GenotypeRefProbabilityPair>& maternal,
                     const std::vector<GenotypeRefProbabilityPair>& paternal)
{
    assert(!maternal.empty() && !paternal.empty());
    return maternal.front().genotype.get().ploidy() == paternal.front().genotype.get().ploidy();
}

auto reduce(std::vector<GenotypeRefProbabilityPair>& likelihoods,
            const PopulationPriorModel& prior_model,
            const TrioModel::Options& options)
{
    auto posteriors = compute_posteriors(likelihoods, prior_model);
    const auto last_posterior_join = reduce(posteriors, get_sample_reduction_count(options.max_joint_genotypes),
                                            options.max_mass_loss);
    assert(last_posterior_join <= std::end(posteriors));
    if (last_posterior_join != std::end(posteriors)) {
        std::sort(std::begin(likelihoods), std::end(likelihoods), GenotypeRefProbabilityPairGenotypeLess {});
        std::sort(std::begin(posteriors), last_posterior_join, GenotypeRefProbabilityPairGenotypeLess {});
        std::vector<GenotypeRefProbabilityPair> reordered_likelihoods {};
        reordered_likelihoods.reserve(posteriors.size());
        std::set_intersection(std::cbegin(likelihoods), std::cend(likelihoods),
                              std::begin(posteriors), last_posterior_join,
                              std::back_inserter(reordered_likelihoods),
                              GenotypeRefProbabilityPairGenotypeLess {});
        std::set_intersection(std::cbegin(likelihoods), std::cend(likelihoods),
                              last_posterior_join, std::end(posteriors),
                              std::back_inserter(reordered_likelihoods),
                              GenotypeRefProbabilityPairGenotypeLess {});
        likelihoods = std::move(reordered_likelihoods);
        auto last_join = std::next(std::cbegin(likelihoods), std::distance(std::begin(posteriors), last_posterior_join));
        return make_reduction_map(likelihoods, last_join, options);
    } else {
        return make_reduction_map(likelihoods, std::cend(likelihoods), options);
    }
}

double probability_of_parents(const Genotype<Haplotype>& mother,
                              const Genotype<Haplotype>& father,
                              const PopulationPriorModel& model)
{
    const std::vector<std::reference_wrapper<const Genotype<Haplotype>>> parental_genotypes {mother, father};
    return model.evaluate(parental_genotypes);
}

template <typename T>
auto size(const ReducedVectorMap<T>& map) noexcept
{
    return static_cast<std::size_t>(std::distance(map.first, map.last));
}

auto join(const ReducedVectorMap<GenotypeRefProbabilityPair>& maternal,
          const ReducedVectorMap<GenotypeRefProbabilityPair>& paternal,
          const PopulationPriorModel& model)
{
    std::vector<ParentsProbabilityPair> result {};
    result.reserve(size(maternal) * size(paternal));
    std::for_each(maternal.first, maternal.last_to_join, [&] (const auto& m) {
        std::for_each(paternal.first, paternal.last_to_join, [&] (const auto& p) {
            result.push_back({m.genotype, p.genotype,
                              m.probability + p.probability + probability_of_parents(m.genotype, p.genotype, model)});
        });
    });
    std::for_each(maternal.last_to_join, maternal.last, [&] (const auto& m) {
        std::for_each(paternal.first, paternal.last_to_partially_join, [&] (const auto& p) {
            result.push_back({m.genotype, p.genotype,
                              m.probability + p.probability + probability_of_parents(m.genotype, p.genotype, model)});
        });
    });
    std::for_each(paternal.last_to_join, paternal.last, [&] (const auto& p) {
        std::for_each(maternal.first, maternal.last_to_partially_join, [&] (const auto& m) {
            result.push_back({m.genotype, p.genotype,
                              m.probability + p.probability + probability_of_parents(m.genotype, p.genotype, model)});
        });
    });
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

auto join(const ReducedVectorMap<ParentsProbabilityPair>& parents,
          const ReducedVectorMap<GenotypeRefProbabilityPair>& child,
          const DeNovoModel& mutation_model)
{
    std::vector<JointProbability> result {};
    result.reserve(size(parents) * size(child));
    std::for_each(parents.first, parents.last_to_join, [&] (const auto& p) {
        std::for_each(child.first, child.last_to_join, [&] (const auto& c) {
            result.push_back({p.maternal, p.paternal, c.genotype, joint_probability(p, c, mutation_model)});
        });
    });
    std::for_each(parents.last_to_join, parents.last, [&] (const auto& p) {
        std::for_each(child.first, child.last_to_partially_join, [&] (const auto& c) {
            result.push_back({p.maternal, p.paternal, c.genotype, joint_probability(p, c, mutation_model)});
        });
    });
    std::for_each(child.last_to_join, child.last, [&] (const auto& c) {
        std::for_each(parents.first, parents.last_to_partially_join, [&] (const auto& p) {
            result.push_back({p.maternal, p.paternal, c.genotype, joint_probability(p, c, mutation_model)});
        });
    });
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
    haplotype_likelihoods.prime(trio_.mother());
    auto maternal_likelihoods = compute_likelihoods(maternal_genotypes, likelihood_model);
    haplotype_likelihoods.prime(trio_.father());
    auto paternal_likelihoods = compute_likelihoods(paternal_genotypes, likelihood_model);
    if (debug_log_) {
        debug::print(stream(*debug_log_), "maternal", maternal_likelihoods);
        debug::print(stream(*debug_log_), "paternal", paternal_likelihoods);
    }
    const auto reduced_maternal_likelihoods = reduce(maternal_likelihoods, prior_model_, options_);
    const auto reduced_paternal_likelihoods = reduce(paternal_likelihoods, prior_model_, options_);
    auto parental_likelihoods = join(reduced_maternal_likelihoods, reduced_paternal_likelihoods, prior_model_);
    clear(maternal_likelihoods);
    clear(paternal_likelihoods);
    if (debug_log_) debug::print(stream(*debug_log_), parental_likelihoods);
    const auto reduced_parental_likelihoods = reduce(parental_likelihoods, options_);
    haplotype_likelihoods.prime(trio_.child());
    auto child_likelihoods = compute_likelihoods(child_genotypes, likelihood_model);
    if (debug_log_) debug::print(stream(*debug_log_), "child", child_likelihoods);
    const auto reduced_child_likelihoods = reduce(child_likelihoods, prior_model_, options_);
    auto joint_likelihoods = join(reduced_parental_likelihoods, reduced_child_likelihoods, mutation_model_);
    clear(parental_likelihoods);
    clear(child_likelihoods);
    const auto evidence = normalise_exp(joint_likelihoods);
    if (debug_log_) debug::print(stream(*debug_log_), joint_likelihoods);
    return {std::move(joint_likelihoods), evidence};
}

TrioModel::InferredLatents
TrioModel::evaluate(const GenotypeVector& genotypes, const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    return evaluate(genotypes, genotypes, genotypes, haplotype_likelihoods);
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
