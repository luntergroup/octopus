// Copyright (c) 2015-2018 Daniel Cooke
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

#include "timers.hpp"

namespace octopus { namespace model {

unsigned TrioModel::max_ploidy() noexcept
{
    return 3;
}

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

const PopulationPriorModel& TrioModel::prior_model() const noexcept
{
    return prior_model_;
}

namespace {

template <typename Container>
void clear(Container& c)
{
    c.clear();
    c.shrink_to_fit();
}

using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;
using GenotypeIndiceVector = GenotypeIndex;
using GenotypeIndiceVectorReference = std::reference_wrapper<const GenotypeIndiceVector>;

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
    const GenotypeIndiceVector* indices = nullptr;
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
    double probability, maternal_likelihood, paternal_likelihood;
    const GenotypeIndiceVector* maternal_indices = nullptr, *paternal_indices = nullptr;
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
    if (prior_model.is_primed()) {
        std::transform(std::cbegin(likelihoods), std::cend(likelihoods), std::begin(posteriors),
                       [&] (const auto& p) {
                           return p.probability + prior_model.evaluate(std::vector<GenotypeIndiceVectorReference> {*p.indices});
                       });
    } else {
        std::transform(std::cbegin(likelihoods), std::cend(likelihoods), std::begin(posteriors),
                       [&] (const auto& p) {
                           return p.probability + prior_model.evaluate(std::vector<GenotypeReference> {p.genotype});
                       });
    }
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
auto compute_posteriors(Iterator first, Iterator last)
{
    std::vector<double> result(std::distance(first, last));
    std::transform(first, last, std::begin(result), [] (const auto& p) { return p.probability; });
    maths::normalise_exp(result);
    return result;
}

template <typename Iterator>
double compute_lost_posterior_mass(Iterator first, Iterator first_removed, Iterator last)
{
    auto posteriors = compute_posteriors(first, last);
    const auto num_removed = std::distance(first, first_removed);
    return std::accumulate(std::next(std::cbegin(posteriors), num_removed), std::cend(posteriors), 0.0);
}

template <typename T>
auto compute_min_joint(const std::vector<T>& zipped, const double max_mass_loss)
{
    assert(!zipped.empty());
    auto posteriors = compute_posteriors(std::cbegin(zipped), std::cend(zipped));
    // Do things in reverse order as floating points have better precision near zero
    std::partial_sum(std::crbegin(posteriors), std::crend(posteriors), std::rbegin(posteriors));
    const auto iter = std::upper_bound(std::crbegin(posteriors), std::crend(posteriors), max_mass_loss);
    const auto num_to_keep = static_cast<std::size_t>(std::distance(iter, std::crend(posteriors)));
    assert(num_to_keep <= zipped.size());
    return std::max(num_to_keep, std::size_t {1});
}

template <typename T>
typename std::vector<T>::iterator
reduce(std::vector<T>& zipped, std::size_t max_joint, const double max_mass_loss, boost::optional<bool&> overflow = boost::none)
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
        if (lost_mass < max_mass_loss) {
            const auto min_joint = compute_min_joint(zipped, max_mass_loss);
            assert(min_joint <= max_joint);
            result = std::next(std::begin(zipped), min_joint);
        } else if (overflow && lost_mass > max_mass_loss) {
            *overflow = true;
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
    auto last_to_join = reduce(zipped, get_sample_reduction_count(options.max_joint_genotypes),
                               options.max_individual_mass_loss);
    return make_reduction_map(zipped, last_to_join, options);
}

struct UniformPriorJointProbabilityHelper
{
    std::size_t index;
    double probability;
};

bool operator==(const UniformPriorJointProbabilityHelper& lhs, const UniformPriorJointProbabilityHelper& rhs) noexcept
{
    return lhs.probability == rhs.probability;
}
bool operator!=(const UniformPriorJointProbabilityHelper& lhs, const UniformPriorJointProbabilityHelper& rhs) noexcept
{
    return lhs.probability != rhs.probability;
}
bool operator<(const UniformPriorJointProbabilityHelper& lhs, const UniformPriorJointProbabilityHelper& rhs) noexcept
{
    return lhs.probability < rhs.probability;
}
bool operator>(const UniformPriorJointProbabilityHelper& lhs, const UniformPriorJointProbabilityHelper& rhs) noexcept
{
    return lhs.probability > rhs.probability;
}

auto reduce(std::vector<ParentsProbabilityPair>& zipped, const TrioModel::Options& options)
{
    const auto reduction_count = get_sample_reduction_count(options.max_joint_genotypes);
    auto last_to_join = reduce(zipped, reduction_count, options.max_joint_mass_loss);
    if (last_to_join != std::cend(zipped)) {
        std::vector<UniformPriorJointProbabilityHelper> likelihood_zipped(zipped.size());
        for (std::size_t i {0}; i < zipped.size(); ++i) {
            likelihood_zipped[i].index = i;
            likelihood_zipped[i].probability = zipped[i].maternal_likelihood + zipped[i].paternal_likelihood;
        }
        likelihood_zipped.erase(reduce(likelihood_zipped, reduction_count, options.max_joint_mass_loss),
                                std::cend(likelihood_zipped));
        std::deque<std::size_t> new_indices {};
        const auto num_posterior_joined = static_cast<std::size_t>(std::distance(std::begin(zipped), last_to_join));
        for (const auto& p : likelihood_zipped) {
            if (p.index >= num_posterior_joined) {
                new_indices.push_back(p.index);
            }
        }
        likelihood_zipped.clear();
        if (!new_indices.empty()) {
            std::sort(std::begin(new_indices), std::end(new_indices));
            auto curr_index = num_posterior_joined;
            // Can use std::partition as it must do a single left-to-right pass due to ForwardIterator
            // and complexity requirements
            last_to_join = std::partition(last_to_join, std::end(zipped), [&] (const auto& p) {
                if (!new_indices.empty() && curr_index++ == new_indices.front()) {
                    new_indices.pop_front();
                    return true;
                } else {
                    return false;
                }
            });
            assert(new_indices.empty());
        }
    }
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
    const auto reduction_count = get_sample_reduction_count(options.max_joint_genotypes);
    bool overflow {false};
    const auto last_likelihood_join = reduce(likelihoods, reduction_count, options.max_individual_mass_loss, overflow);
    if (last_likelihood_join == std::end(likelihoods)) {
        return make_reduction_map(likelihoods, std::cend(likelihoods), options);
    } else if (!overflow) {
        return make_reduction_map(likelihoods, last_likelihood_join, options);
    } else {
        auto posteriors = compute_posteriors(likelihoods, prior_model);
        const auto last_posterior_join = reduce(posteriors, reduction_count, options.max_individual_mass_loss);
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
}

double joint_probability(const Genotype<Haplotype>& mother, const Genotype<Haplotype>& father,
                         const PopulationPriorModel& model)
{
    const std::vector<std::reference_wrapper<const Genotype<Haplotype>>> parental_genotypes {mother, father};
    return model.evaluate(parental_genotypes);
}

double joint_probability(const GenotypeIndiceVector& mother, const GenotypeIndiceVector& father,
                         const PopulationPriorModel& model)
{
    const std::vector<GenotypeIndiceVectorReference> parental_genotypes {mother, father};
    return model.evaluate(parental_genotypes);
}

double joint_probability(const GenotypeRefProbabilityPair& mother, const GenotypeRefProbabilityPair& father,
                         const PopulationPriorModel& model)
{
    if (mother.indices && father.indices) {
        return mother.probability + father.probability + joint_probability(*mother.indices, *father.indices, model);
    } else {
        return mother.probability + father.probability + joint_probability(mother.genotype, father.genotype, model);
    }
}

template <typename T1, typename T2>
auto join_size(const ReducedVectorMap<T1>& first, const ReducedVectorMap<T2>& second) noexcept
{
    using std::distance;
    std::size_t result {0};
    result += distance(first.first, first.last_to_join) * distance(second.first, second.last_to_join);
    result += distance(first.last_to_join, first.last) * distance(second.first, second.last_to_partially_join);
    result += distance(second.last_to_join, second.last) * distance(first.first, first.last_to_partially_join);
    return result;
}

auto join(const ReducedVectorMap<GenotypeRefProbabilityPair>& maternal,
          const ReducedVectorMap<GenotypeRefProbabilityPair>& paternal,
          const PopulationPriorModel& model)
{
    std::vector<ParentsProbabilityPair> result {};
    result.reserve(join_size(maternal, paternal));
    std::for_each(maternal.first, maternal.last_to_join, [&] (const auto& m) {
        std::for_each(paternal.first, paternal.last_to_join, [&] (const auto& p) {
            result.push_back({m.genotype, p.genotype, joint_probability(m, p, model),
                              m.probability, m.probability, m.indices, p.indices});
        });
    });
    std::for_each(maternal.last_to_join, maternal.last, [&] (const auto& m) {
        std::for_each(paternal.first, paternal.last_to_partially_join, [&] (const auto& p) {
            result.push_back({m.genotype, p.genotype, joint_probability(m, p, model),
                              m.probability, m.probability, m.indices, p.indices});
        });
    });
    std::for_each(paternal.last_to_join, paternal.last, [&] (const auto& p) {
        std::for_each(maternal.first, maternal.last_to_partially_join, [&] (const auto& m) {
            result.push_back({m.genotype, p.genotype, joint_probability(m, p, model),
                              m.probability, m.probability, m.indices, p.indices});
        });
    });
    return result;
}

bool is_haploid(const GenotypeIndex& genotype) noexcept
{
    return genotype.size() == 1;
}

bool is_diploid(const GenotypeIndex& genotype) noexcept
{
    return genotype.size() == 2;
}

bool is_triploid(const GenotypeIndex& genotype) noexcept
{
    return genotype.size() == 3;
}

template <typename G>
bool all_haploid(const G& child, const G& mother, const G& father) noexcept
{
    return is_haploid(child) && is_haploid(mother) && is_haploid(father);
}

template <typename G>
bool all_diploid(const G& child, const G& mother, const G& father) noexcept
{
    return is_diploid(child) && is_diploid(mother) && is_diploid(father);
}

template <typename G>
bool all_triploid(const G& child, const G& mother, const G& father) noexcept
{
    return is_triploid(child) && is_triploid(mother) && is_triploid(father);
}

template <typename H, typename G>
double probability_of_child_given_haploid_parent(const H& child, const G& parent,
                                                 const DeNovoModel& mutation_model)
{
    return mutation_model.evaluate(child, parent[0]);
}

template <typename H, typename G>
double probability_of_child_given_diploid_parent(const H& child, const G& parent,
                                                 const DeNovoModel& mutation_model)
{
    static const double ln2 {std::log(2)};
    const auto p1 = mutation_model.evaluate(child, parent[0]);
    const auto p2 = mutation_model.evaluate(child, parent[1]);
    return maths::log_sum_exp(p1, p2) - ln2;
}

template <typename H, typename G>
double probability_of_child_given_diploid_parents(const H& child_from_mother, const H& child_from_father,
                                                  const G& mother, const G& father,
                                                  const DeNovoModel& mutation_model)
{
    return probability_of_child_given_diploid_parent(child_from_mother, mother, mutation_model)
            + probability_of_child_given_diploid_parent(child_from_father, father, mutation_model);
}

template <typename H, typename G>
double probability_of_child_given_triploid_parent(const H& child, const G& parent,
                                                  const DeNovoModel& mutation_model)
{
    static const double ln3 {std::log(3)};
    const auto p1 = mutation_model.evaluate(child, parent[0]);
    const auto p2 = mutation_model.evaluate(child, parent[1]);
    const auto p3 = mutation_model.evaluate(child, parent[2]);
    return maths::log_sum_exp(p1, p2, p3) - ln3;
}

template <typename H, typename G>
double probability_of_child_given_triploid_parents(const H& child1_from_mother,
                                                   const H& child2_from_mother,
                                                   const H& child_from_father,
                                                   const G& mother, const G& father,
                                                   const DeNovoModel& mutation_model)
{
    return probability_of_child_given_triploid_parent(child1_from_mother, mother, mutation_model)
           + probability_of_child_given_triploid_parent(child2_from_mother, mother, mutation_model)
           + probability_of_child_given_triploid_parent(child_from_father, father, mutation_model);
}

template <unsigned ChildPloidy, unsigned MotherPloidy, unsigned FatherPloidy>
struct ProbabilityOfChildGivenParents
{
    ProbabilityOfChildGivenParents(const DeNovoModel& mutation_model) : mutation_model {mutation_model} {}
    
    template <typename G>
    double operator()(const G& child, const G& mother, const G& father)
    {
        return 0;
    }
    
    const DeNovoModel& mutation_model;
};

template <> struct ProbabilityOfChildGivenParents<2, 2, 2>
{
    ProbabilityOfChildGivenParents(const DeNovoModel& mutation_model) : mutation_model {mutation_model} {}
    
    template <typename G>
    double operator()(const G& child, const G& mother, const G& father)
    {
        static const double ln2 {std::log(2)};
        const auto p1 = probability_of_child_given_diploid_parents(child[0], child[1], mother, father, mutation_model);
        const auto p2 = probability_of_child_given_diploid_parents(child[1], child[0], mother, father, mutation_model);
        return maths::log_sum_exp(p1, p2) - ln2;
    }
    
    const DeNovoModel& mutation_model;
};

template <> struct ProbabilityOfChildGivenParents<3, 3, 3>
{
    ProbabilityOfChildGivenParents(const DeNovoModel& mutation_model) : mutation_model {mutation_model} {}
    
    template <typename G>
    double operator()(const G& child, const G& mother, const G& father)
    {
        static const double ln6 {std::log(6)};
        const auto p1 = probability_of_child_given_triploid_parents(child[0], child[1], child[2], mother, father, mutation_model);
        const auto p2 = probability_of_child_given_triploid_parents(child[0], child[2], child[1], mother, father, mutation_model);
        const auto p3 = probability_of_child_given_triploid_parents(child[1], child[0], child[2], mother, father, mutation_model);
        const auto p4 = probability_of_child_given_triploid_parents(child[1], child[2], child[0], mother, father, mutation_model);
        const auto p5 = probability_of_child_given_triploid_parents(child[2], child[0], child[1], mother, father, mutation_model);
        const auto p6 = probability_of_child_given_triploid_parents(child[2], child[1], child[0], mother, father, mutation_model);
        return maths::log_sum_exp({p1, p2, p3, p4, p5, p6}) - ln6;
    }
    
    const DeNovoModel& mutation_model;
};

template <> struct ProbabilityOfChildGivenParents<2, 2, 1>
{
    ProbabilityOfChildGivenParents(const DeNovoModel& mutation_model) : mutation_model {mutation_model} {}
    
    template <typename G>
    double operator()(const G& child, const G& mother, const G& father)
    {
        static const double ln2 {std::log(2)};
        const auto p1 = probability_of_child_given_diploid_parent(child[0], mother, mutation_model);
        const auto p2 = probability_of_child_given_diploid_parent(child[1], mother, mutation_model);
        const auto p3 = probability_of_child_given_haploid_parent(child[0], father, mutation_model);
        const auto p4 = probability_of_child_given_haploid_parent(child[1], father, mutation_model);
        return maths::log_sum_exp(p1 + p4, p2 + p3) - ln2;
    }
    
    const DeNovoModel& mutation_model;
};

template <> struct ProbabilityOfChildGivenParents<1, 2, 1>
{
    ProbabilityOfChildGivenParents(const DeNovoModel& mutation_model) : mutation_model {mutation_model} {}
    
    template <typename G>
    double operator()(const G& child, const G& mother, const G& father)
    {
        return probability_of_child_given_haploid_parent(child[0], father, mutation_model);
    }
    
    const DeNovoModel& mutation_model;
};

using JointProbability = TrioModel::Latents::JointProbability;

template <typename F>
auto joint_probability(const ParentsProbabilityPair& parents,
                       const GenotypeRefProbabilityPair& child,
                       F joint_probability_function)
{
    if (child.indices && parents.maternal_indices && parents.paternal_indices) {
        return parents.probability + child.probability
               + joint_probability_function(*child.indices, *parents.maternal_indices, *parents.paternal_indices);
    } else {
        return parents.probability + child.probability
               + joint_probability_function(child.genotype.get(), parents.maternal.get(), parents.paternal.get());
    }
}

template <typename F>
auto join(const ReducedVectorMap<ParentsProbabilityPair>& parents,
          const ReducedVectorMap<GenotypeRefProbabilityPair>& child,
          F jpdf)
{
    std::vector<JointProbability> result {};
    result.reserve(join_size(parents, child));
    std::for_each(parents.first, parents.last_to_join, [&] (const auto& p) {
        std::for_each(child.first, child.last_to_join, [&] (const auto& c) {
            result.push_back({p.maternal, p.paternal, c.genotype, joint_probability(p, c, jpdf)});
        });
    });
    std::for_each(parents.last_to_join, parents.last, [&] (const auto& p) {
        std::for_each(child.first, child.last_to_partially_join, [&] (const auto& c) {
            result.push_back({p.maternal, p.paternal, c.genotype, joint_probability(p, c, jpdf)});
        });
    });
    std::for_each(child.last_to_join, child.last, [&] (const auto& c) {
        std::for_each(parents.first, parents.last_to_partially_join, [&] (const auto& p) {
            result.push_back({p.maternal, p.paternal, c.genotype, joint_probability(p, c, jpdf)});
        });
    });
    return result;
}

auto join(const ReducedVectorMap<ParentsProbabilityPair>& parents,
          const ReducedVectorMap<GenotypeRefProbabilityPair>& child,
          const DeNovoModel& mutation_model)
{
    const auto maternal_ploidy = parents.first->maternal.get().ploidy();
    const auto paternal_ploidy = parents.first->paternal.get().ploidy();
    const auto child_ploidy    = child.first->genotype.get().ploidy();
    if (child_ploidy == 1) {
        if (paternal_ploidy == 1) {
            return join(parents, child, ProbabilityOfChildGivenParents<1, 2, 1> {mutation_model});
        }
    } else if (child_ploidy == 2) {
        if (maternal_ploidy == 2) {
            if (paternal_ploidy == 1) {
                return join(parents, child, ProbabilityOfChildGivenParents<2, 2, 1> {mutation_model});
            }
            if (paternal_ploidy == 2) {
                return join(parents, child, ProbabilityOfChildGivenParents<2, 2, 2> {mutation_model});
            }
        } else {
        
        }
    } else if (child_ploidy == 3 && maternal_ploidy == 3 && paternal_ploidy == 3) {
        return join(parents, child, ProbabilityOfChildGivenParents<3, 3, 3> {mutation_model});
    }
    throw std::runtime_error {"TrioModel: unimplemented joint probability function"};
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
    if (maternal_genotypes.empty()) {
        haplotype_likelihoods.prime(trio_.father());
        return evaluate_allosome(paternal_genotypes, child_genotypes, haplotype_likelihoods);
    }
    if (paternal_genotypes.empty()) {
        haplotype_likelihoods.prime(trio_.mother());
        return evaluate_allosome(maternal_genotypes, child_genotypes, haplotype_likelihoods);
    }
    assert(!maternal_genotypes.empty() && !paternal_genotypes.empty() && !child_genotypes.empty());
    const GermlineLikelihoodModel likelihood_model {haplotype_likelihoods};
    haplotype_likelihoods.prime(trio_.mother());
    auto maternal_likelihoods = compute_likelihoods(maternal_genotypes, likelihood_model);
    haplotype_likelihoods.prime(trio_.father());
    auto paternal_likelihoods = compute_likelihoods(paternal_genotypes, likelihood_model);
    haplotype_likelihoods.prime(trio_.child());
    auto child_likelihoods = compute_likelihoods(child_genotypes, likelihood_model);
    if (debug_log_) {
        debug::print(stream(*debug_log_), "maternal", maternal_likelihoods);
        debug::print(stream(*debug_log_), "paternal", paternal_likelihoods);
        debug::print(stream(*debug_log_), "child", child_likelihoods);
    }
    const auto reduced_maternal_likelihoods = reduce(maternal_likelihoods, prior_model_, options_);
    const auto reduced_paternal_likelihoods = reduce(paternal_likelihoods, prior_model_, options_);
    const auto reduced_child_likelihoods    = reduce(child_likelihoods, prior_model_, options_);
    auto parental_likelihoods = join(reduced_maternal_likelihoods, reduced_paternal_likelihoods, prior_model_);
    if (debug_log_) debug::print(stream(*debug_log_), parental_likelihoods);
    const auto reduced_parental_likelihoods = reduce(parental_likelihoods, options_);
    auto joint_likelihoods = join(reduced_parental_likelihoods, reduced_child_likelihoods, mutation_model_);
    if (debug_log_) debug::print(stream(*debug_log_), joint_likelihoods);
    const auto evidence = normalise_exp(joint_likelihoods);
    return {std::move(joint_likelihoods), evidence};
}

TrioModel::InferredLatents
TrioModel::evaluate(const GenotypeVector& genotypes, const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    return evaluate(genotypes, genotypes, genotypes, haplotype_likelihoods);
}

TrioModel::InferredLatents
TrioModel::evaluate(const GenotypeVector& genotypes, std::vector<GenotypeIndex>& genotype_indices,
                    const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(prior_model_.is_primed() && mutation_model_.is_primed());
    const GermlineLikelihoodModel likelihood_model {haplotype_likelihoods};
    haplotype_likelihoods.prime(trio_.mother());
    auto maternal_likelihoods = compute_likelihoods(genotypes, likelihood_model);
    haplotype_likelihoods.prime(trio_.father());
    auto paternal_likelihoods = compute_likelihoods(genotypes, likelihood_model);
    haplotype_likelihoods.prime(trio_.child());
    auto child_likelihoods = compute_likelihoods(genotypes, likelihood_model);
    for (std::size_t i {0}; i < genotypes.size(); ++i) {
        maternal_likelihoods[i].indices = &genotype_indices[i];
        paternal_likelihoods[i].indices = &genotype_indices[i];
        child_likelihoods[i].indices    = &genotype_indices[i];
    }
    if (debug_log_) {
        debug::print(stream(*debug_log_), "maternal", maternal_likelihoods);
        debug::print(stream(*debug_log_), "paternal", paternal_likelihoods);
        debug::print(stream(*debug_log_), "child", child_likelihoods);
    }
    const auto reduced_maternal_likelihoods = reduce(maternal_likelihoods, prior_model_, options_);
    const auto reduced_paternal_likelihoods = reduce(paternal_likelihoods, prior_model_, options_);
    const auto reduced_child_likelihoods    = reduce(child_likelihoods, prior_model_, options_);
    auto parental_likelihoods = join(reduced_maternal_likelihoods, reduced_paternal_likelihoods, prior_model_);
    if (debug_log_) debug::print(stream(*debug_log_), parental_likelihoods);
    const auto reduced_parental_likelihoods = reduce(parental_likelihoods, options_);
    auto joint_likelihoods = join(reduced_parental_likelihoods, reduced_child_likelihoods, mutation_model_);
    if (debug_log_) debug::print(stream(*debug_log_), joint_likelihoods);
    const auto evidence = normalise_exp(joint_likelihoods);
    return {std::move(joint_likelihoods), evidence};
}

double probability_of_child_given_parent(const Genotype<Haplotype>& child,
                                         const Genotype<Haplotype>& parent,
                                         const DeNovoModel& mutation_model)
{
    if (is_haploid(child) && is_haploid(parent)) {
        return mutation_model.evaluate(child[0], parent[0]);
    }
    return 0; // TODO
}

auto joint_probability(const GenotypeRefProbabilityPair& parent,
                       const GenotypeRefProbabilityPair& child,
                       const DeNovoModel& mutation_model)
{
    return parent.probability + child.probability
           + probability_of_child_given_parent(child.genotype, parent.genotype, mutation_model);
}

auto join(const ReducedVectorMap<GenotypeRefProbabilityPair>& parent,
          const ReducedVectorMap<GenotypeRefProbabilityPair>& child,
          const DeNovoModel& mutation_model)
{
    std::vector<JointProbability> result {};
    result.reserve(join_size(parent, child));
    std::for_each(parent.first, parent.last_to_join, [&] (const auto& p) {
        std::for_each(child.first, child.last_to_join, [&] (const auto& c) {
            result.push_back({p.genotype, p.genotype, c.genotype, joint_probability(p, c, mutation_model)});
        });
    });
    std::for_each(parent.last_to_join, parent.last, [&] (const auto& p) {
        std::for_each(child.first, child.last_to_partially_join, [&] (const auto& c) {
            result.push_back({p.genotype, p.genotype, c.genotype, joint_probability(p, c, mutation_model)});
        });
    });
    std::for_each(child.last_to_join, child.last, [&] (const auto& c) {
        std::for_each(parent.first, parent.last_to_partially_join, [&] (const auto& p) {
            result.push_back({p.genotype, p.genotype, c.genotype, joint_probability(p, c, mutation_model)});
        });
    });
    return result;
}

TrioModel::InferredLatents
TrioModel::evaluate_allosome(const GenotypeVector& parent_genotypes,
                             const GenotypeVector& child_genotypes,
                             const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!parent_genotypes.empty() && child_genotypes.empty());
    const GermlineLikelihoodModel likelihood_model {haplotype_likelihoods};
    assert(haplotype_likelihoods.is_primed());
    auto parent_likelihoods = compute_likelihoods(parent_genotypes, likelihood_model);
    if (debug_log_) debug::print(stream(*debug_log_), "parent", parent_likelihoods);
    const auto reduced_parent_likelihoods = reduce(parent_likelihoods, prior_model_, options_);
    haplotype_likelihoods.prime(trio_.child());
    auto child_likelihoods = compute_likelihoods(child_genotypes, likelihood_model);
    if (debug_log_) debug::print(stream(*debug_log_), "child", child_likelihoods);
    const auto reduced_child_likelihoods = reduce(child_likelihoods, prior_model_, options_);
    auto joint_likelihoods = join(reduced_parent_likelihoods, reduced_child_likelihoods, mutation_model_);
    clear(parent_likelihoods);
    clear(child_likelihoods);
    const auto evidence = normalise_exp(joint_likelihoods);
    if (debug_log_) debug::print(stream(*debug_log_), joint_likelihoods);
    return {std::move(joint_likelihoods), evidence};
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
