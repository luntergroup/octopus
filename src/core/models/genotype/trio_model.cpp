// Copyright (c) 2015-2019 Daniel Cooke
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

#include <boost/iterator/transform_iterator.hpp>

#include "utils/maths.hpp"
#include "constant_mixture_genotype_likelihood_model.hpp"

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
                         const ConstantMixtureGenotypeLikelihoodModel& model)
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

struct ProbabilityGetter
{
    template <typename T> const auto& operator()(const T& x) const noexcept { return x.probability; }
};

template <typename Iterator>
auto compute_lost_posterior_mass(Iterator first, Iterator first_removed, Iterator last)
{
    using boost::make_transform_iterator;
    const auto first_prob = make_transform_iterator(first, ProbabilityGetter {});
    const auto first_removed_prob = make_transform_iterator(first_removed, ProbabilityGetter {});
    const auto last_prob = make_transform_iterator(last, ProbabilityGetter {});
    return maths::log_sum_exp(first_removed_prob, last_prob) - maths::log_sum_exp(first_prob, last_prob);
}

template <typename T>
auto compute_min_joint(const std::vector<T>& zipped, const double max_log_mass_loss, boost::optional<double>& lost_log_mass)
{
    assert(!zipped.empty());
    using boost::make_transform_iterator;
    std::vector<double> cum_log_probs(make_transform_iterator(std::cbegin(zipped), ProbabilityGetter {}),
                                      make_transform_iterator(std::cend(zipped), ProbabilityGetter {}));
    maths::normalise_logs(cum_log_probs);
    const static auto log_sum_exp = [] (auto sum, auto x) { return maths::log_sum_exp(sum, x); };
    std::partial_sum(std::crbegin(cum_log_probs), std::crend(cum_log_probs), std::rbegin(cum_log_probs), log_sum_exp);
    const auto min_itr = std::upper_bound(std::crbegin(cum_log_probs), std::crend(cum_log_probs), max_log_mass_loss);
    const auto num_to_keep = static_cast<std::size_t>(std::distance(min_itr, std::crend(cum_log_probs)));
    assert(num_to_keep <= zipped.size());
    if (num_to_keep < zipped.size()) {
        assert(min_itr != std::crbegin(cum_log_probs));
        lost_log_mass = *std::prev(min_itr);
    }
    return std::max(num_to_keep, std::size_t {1});
}

template <typename T>
typename std::vector<T>::iterator
reduce(std::vector<T>& zipped, std::size_t max_joint,
       const double max_log_probability_loss,
       boost::optional<double>& lost_log_mass)
{
    max_joint = std::min(max_joint, zipped.size());
    if (all_equal(zipped)) {
        static std::default_random_engine gen {};
        std::shuffle(std::begin(zipped), std::end(zipped), gen);
        const auto min_joint = static_cast<std::size_t>(std::floor(zipped.size() * (1 - std::exp(max_log_probability_loss))));
        lost_log_mass = boost::none;
        return std::next(std::begin(zipped), std::min(max_joint, min_joint));
    } else {
        auto result = std::next(std::begin(zipped), max_joint);
        std::partial_sort(std::begin(zipped), result, std::end(zipped), std::greater<> {});
        if (max_joint < zipped.size()) {
            lost_log_mass = compute_lost_posterior_mass(std::begin(zipped), result, std::end(zipped));
        }
        if (!lost_log_mass || *lost_log_mass < max_log_probability_loss) {
            auto num_joint = compute_min_joint(zipped, max_log_probability_loss, lost_log_mass);
            assert(num_joint <= max_joint);
            result = std::next(std::begin(zipped), num_joint);
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
    auto last_to_join = std::cend(zipped);
    if (options.max_genotype_combinations) {
        last_to_join = reduce(zipped, get_sample_reduction_count(*options.max_genotype_combinations), options.max_individual_log_probability_loss);
    }
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

template <typename Iterator, typename UnaryOperation>
auto index_zip(Iterator first, const Iterator last, UnaryOperation op)
{
    using T = typename std::iterator_traits<Iterator>::value_type;
    using R = std::result_of_t<UnaryOperation(const T&)>;
    const auto n = static_cast<std::size_t>(std::distance(first, last));
    std::vector<std::pair<R, std::size_t>> result(n);
    for (std::size_t idx {0}; idx < n; ++idx) {
        result[idx] = std::make_pair(op(*first++), idx);
    }
    return result;
}

template <typename Iterator>
auto order_indices_by_probability(const Iterator first, const Iterator last, std::size_t k)
{
    auto zipped = index_zip(first, last, [] (const auto& pair) { return pair.probability; });
    k = std::min(zipped.size(), k);
    const auto kth = std::next(std::begin(zipped), k);
    std::partial_sort(std::begin(zipped), kth, std::end(zipped), std::greater<> {});
    std::vector<std::size_t> result(k);
    std::transform(std::begin(zipped), kth, std::begin(result), [] (const auto& p) { return p.second; });
    return result;
}

template <typename Iterator>
auto select_top_k_haplotypes(const Iterator first, Iterator last, const std::size_t k)
{
    const auto top_pair_indices = order_indices_by_probability(first, last, k);
    std::set<Haplotype> result {};
    for (const auto idx : top_pair_indices) {
        if (result.size() >= k) break;
        const auto& genotype = std::next(first, idx)->genotype.get();
        std::copy(std::cbegin(genotype), std::cend(genotype), std::inserter(result, std::begin(result)));
    }
    return result;
}

using ChildReductionMap = ReducedVectorMap<GenotypeRefProbabilityPair>;

auto select_top_k_haplotypes(const ChildReductionMap& child, const std::size_t k)
{
    return select_top_k_haplotypes(child.first, child.last_to_join, k);
}

using ParentsProbabilityPairIterator = std::vector<ParentsProbabilityPair>::iterator;

bool is_represented(const Haplotype& haplotype, const ParentsProbabilityPair& parent)
{
    return parent.maternal.get().contains(haplotype) || parent.paternal.get().contains(haplotype);
}

bool is_represented(const Haplotype& haplotype,
                    const ParentsProbabilityPairIterator first_parent,
                    const ParentsProbabilityPairIterator last_parent)
{
    return std::any_of(first_parent, last_parent, [&] (const auto& parent) { return is_represented(haplotype, parent); });
}

auto find_represented(const Haplotype& haplotype,
                      const ParentsProbabilityPairIterator first_parent,
                      const ParentsProbabilityPairIterator last_parent)
{
    return std::find_if(first_parent, last_parent, [&] (const auto& parent) { return is_represented(haplotype, parent); });
}

auto reduce(std::vector<ParentsProbabilityPair>& parents, const ChildReductionMap& child,
            boost::optional<double>& lost_log_mass, const TrioModel::Options& options)
{
    if (!options.max_genotype_combinations) return make_reduction_map(parents, std::cend(parents), options);
    const auto reduction_count = get_sample_reduction_count(*options.max_genotype_combinations);
    auto last_to_join = reduce(parents, reduction_count, options.max_joint_log_probability_loss, lost_log_mass);
    if (last_to_join != std::cend(parents)) {
        std::vector<UniformPriorJointProbabilityHelper> uniform_parents(parents.size());
        for (std::size_t i {0}; i < parents.size(); ++i) {
            uniform_parents[i].index = i;
            uniform_parents[i].probability = parents[i].maternal_likelihood + parents[i].paternal_likelihood;
        }
        uniform_parents.erase(reduce(uniform_parents, reduction_count, options.max_joint_log_probability_loss, lost_log_mass), std::cend(uniform_parents));
        std::deque<std::size_t> new_indices {};
        const auto num_posterior_joined = static_cast<std::size_t>(std::distance(std::begin(parents), last_to_join));
        for (const auto& p : uniform_parents) {
            if (p.index >= num_posterior_joined) {
                new_indices.push_back(p.index);
            }
        }
        uniform_parents.clear();
        if (!new_indices.empty()) {
            std::sort(std::begin(new_indices), std::end(new_indices));
            auto curr_index = num_posterior_joined;
            // Can use std::partition as it must do a single left-to-right pass due to ForwardIterator and complexity requirements
            last_to_join = std::partition(last_to_join, std::end(parents), [&] (const auto& p) {
                if (!new_indices.empty() && curr_index++ == new_indices.front()) {
                    new_indices.pop_front();
                    return true;
                } else {
                    return false;
                }
            });
            assert(new_indices.empty());
        }
        // We want to make sure haplotypes that the child may have are survive to avoid false positive de novo child haplotypes
        const auto top_child_haplotypes = select_top_k_haplotypes(child, 4);
        for (const auto& haplotype : top_child_haplotypes) {
            if (!is_represented(haplotype, std::begin(parents), last_to_join)) {
                const auto first_represented_parent = find_represented(haplotype, last_to_join, std::end(parents));
                if (first_represented_parent != std::end(parents)) {
                    std::iter_swap(first_represented_parent, last_to_join);
                    ++last_to_join;
                }
            }
        }
    }
    return make_reduction_map(parents, last_to_join, options);
}

bool are_same_ploidy(const std::vector<GenotypeRefProbabilityPair>& maternal,
                     const std::vector<GenotypeRefProbabilityPair>& paternal)
{
    assert(!maternal.empty() && !paternal.empty());
    return maternal.front().genotype.get().ploidy() == paternal.front().genotype.get().ploidy();
}

auto reduce(std::vector<GenotypeRefProbabilityPair>& likelihoods,
            const PopulationPriorModel& prior_model,
            boost::optional<double>& lost_log_mass,
            const TrioModel::Options& options)
{
    if (!options.max_genotype_combinations) return make_reduction_map(likelihoods, std::cend(likelihoods), options);
    const auto reduction_count = get_sample_reduction_count(*options.max_genotype_combinations);
    const auto last_likelihood_join = reduce(likelihoods, reduction_count, options.max_individual_log_probability_loss, lost_log_mass);
    if (last_likelihood_join == std::end(likelihoods)) {
        return make_reduction_map(likelihoods, std::cend(likelihoods), options);
    } else if (lost_log_mass && *lost_log_mass > options.max_individual_log_probability_loss) {
        return make_reduction_map(likelihoods, last_likelihood_join, options); // overflow
    } else {
        auto posteriors = compute_posteriors(likelihoods, prior_model);
        const auto last_posterior_join = reduce(posteriors, reduction_count, options.max_individual_log_probability_loss, lost_log_mass);
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
        return probability_of_child_given_diploid_parent(child[0], mother, mutation_model);
    }
    
    const DeNovoModel& mutation_model;
};

template <> struct ProbabilityOfChildGivenParents<1, 0, 1>
{
    ProbabilityOfChildGivenParents(const DeNovoModel& mutation_model) : mutation_model {mutation_model} {}
    
    template <typename G>
    double operator()(const G& child, const G& mother, const G& father)
    {
        return probability_of_child_given_haploid_parent(child[0], father, mutation_model);
    }
    
    const DeNovoModel& mutation_model;
};

template <> struct ProbabilityOfChildGivenParents<1, 1, 1>
{
    ProbabilityOfChildGivenParents(const DeNovoModel& mutation_model) : mutation_model {mutation_model} {}
    
    template <typename G>
    double operator()(const G& child, const G& mother, const G& father)
    {
        static const double ln2 {std::log(2)};
        const auto p1 = probability_of_child_given_haploid_parent(child[0], mother, mutation_model);
        const auto p2 = probability_of_child_given_haploid_parent(child[0], father, mutation_model);
        return maths::log_sum_exp(p1, p2) - ln2;
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
            result.push_back({p.maternal, p.paternal, c.genotype, joint_probability(p, c, jpdf), 0.0});
        });
    });
    std::for_each(parents.last_to_join, parents.last, [&] (const auto& p) {
        std::for_each(child.first, child.last_to_partially_join, [&] (const auto& c) {
            result.push_back({p.maternal, p.paternal, c.genotype, joint_probability(p, c, jpdf), 0.0});
        });
    });
    std::for_each(child.last_to_join, child.last, [&] (const auto& c) {
        std::for_each(parents.first, parents.last_to_partially_join, [&] (const auto& p) {
            result.push_back({p.maternal, p.paternal, c.genotype, joint_probability(p, c, jpdf), 0.0});
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
            if (maternal_ploidy == 0) {
                return join(parents, child, ProbabilityOfChildGivenParents<1, 0, 1> {mutation_model});
            }
            if (maternal_ploidy == 1) {
                return join(parents, child, ProbabilityOfChildGivenParents<1, 1, 1> {mutation_model});
            }
            if (maternal_ploidy == 2) {
                return join(parents, child, ProbabilityOfChildGivenParents<1, 2, 1> {mutation_model});
            }
        }
    } else if (child_ploidy == 2) {
        if (maternal_ploidy == 2) {
            if (paternal_ploidy == 1) {
                return join(parents, child, ProbabilityOfChildGivenParents<2, 2, 1> {mutation_model});
            }
            if (paternal_ploidy == 2) {
                return join(parents, child, ProbabilityOfChildGivenParents<2, 2, 2> {mutation_model});
            }
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
                   std::begin(result), [] (const auto& p) { return p.log_probability; });
    return result;
}

auto normalise_exp(std::vector<JointProbability>& joint_likelihoods)
{
    auto log_likelihoods = extract_probabilities(joint_likelihoods);
    const auto norm = maths::normalise_logs(log_likelihoods);
    auto iter = std::cbegin(log_likelihoods);
    for (auto& p : joint_likelihoods) {
        p.log_probability = *iter++;
        p.probability = std::exp(p.log_probability);
    }
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
                    const HaplotypeLikelihoodArray& haplotype_likelihoods) const
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
    const ConstantMixtureGenotypeLikelihoodModel likelihood_model {haplotype_likelihoods};
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
    boost::optional<double> lost_log_mass {};
    const auto reduced_maternal_likelihoods = reduce(maternal_likelihoods, prior_model_, lost_log_mass, options_);
    const auto reduced_paternal_likelihoods = reduce(paternal_likelihoods, prior_model_, lost_log_mass, options_);
    const auto reduced_child_likelihoods    = reduce(child_likelihoods, prior_model_, lost_log_mass, options_);
    auto parental_likelihoods = join(reduced_maternal_likelihoods, reduced_paternal_likelihoods, prior_model_);
    if (debug_log_) debug::print(stream(*debug_log_), parental_likelihoods);
    const auto reduced_parental_likelihoods = reduce(parental_likelihoods, reduced_child_likelihoods, lost_log_mass, options_);
    auto joint_likelihoods = join(reduced_parental_likelihoods, reduced_child_likelihoods, mutation_model_);
    if (debug_log_) debug::print(stream(*debug_log_), joint_likelihoods);
    const auto evidence = normalise_exp(joint_likelihoods);
    if (lost_log_mass) *lost_log_mass *= 2 * std::distance(reduced_child_likelihoods.first, reduced_child_likelihoods.last_to_join);
    return {std::move(joint_likelihoods), evidence, lost_log_mass};
}

TrioModel::InferredLatents
TrioModel::evaluate(const GenotypeVector& genotypes, const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    return evaluate(genotypes, genotypes, genotypes, haplotype_likelihoods);
}

TrioModel::InferredLatents
TrioModel::evaluate(const GenotypeVector& genotypes,
                    std::vector<GenotypeIndex>& genotype_indices,
                    const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    assert(prior_model_.is_primed() && mutation_model_.is_primed());
    const ConstantMixtureGenotypeLikelihoodModel likelihood_model {haplotype_likelihoods};
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
    boost::optional<double> lost_log_mass {};
    const auto reduced_maternal_likelihoods = reduce(maternal_likelihoods, prior_model_, lost_log_mass, options_);
    const auto reduced_paternal_likelihoods = reduce(paternal_likelihoods, prior_model_, lost_log_mass, options_);
    const auto reduced_child_likelihoods    = reduce(child_likelihoods, prior_model_, lost_log_mass, options_);
    auto parental_likelihoods = join(reduced_maternal_likelihoods, reduced_paternal_likelihoods, prior_model_);
    if (debug_log_) debug::print(stream(*debug_log_), parental_likelihoods);
    const auto reduced_parental_likelihoods = reduce(parental_likelihoods, reduced_child_likelihoods, lost_log_mass, options_);
    auto joint_likelihoods = join(reduced_parental_likelihoods, reduced_child_likelihoods, mutation_model_);
    if (debug_log_) debug::print(stream(*debug_log_), joint_likelihoods);
    const auto evidence = normalise_exp(joint_likelihoods);
    if (lost_log_mass) *lost_log_mass *= 2 * std::distance(reduced_child_likelihoods.first, reduced_child_likelihoods.last_to_join);
    return {std::move(joint_likelihoods), evidence, lost_log_mass};
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
            result.push_back({p.genotype, p.genotype, c.genotype, joint_probability(p, c, mutation_model), 0.0});
        });
    });
    std::for_each(parent.last_to_join, parent.last, [&] (const auto& p) {
        std::for_each(child.first, child.last_to_partially_join, [&] (const auto& c) {
            result.push_back({p.genotype, p.genotype, c.genotype, joint_probability(p, c, mutation_model), 0.0});
        });
    });
    std::for_each(child.last_to_join, child.last, [&] (const auto& c) {
        std::for_each(parent.first, parent.last_to_partially_join, [&] (const auto& p) {
            result.push_back({p.genotype, p.genotype, c.genotype, joint_probability(p, c, mutation_model), 0.0});
        });
    });
    return result;
}

TrioModel::InferredLatents
TrioModel::evaluate_allosome(const GenotypeVector& parent_genotypes,
                             const GenotypeVector& child_genotypes,
                             const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    assert(!parent_genotypes.empty() && !child_genotypes.empty());
    const ConstantMixtureGenotypeLikelihoodModel likelihood_model {haplotype_likelihoods};
    assert(haplotype_likelihoods.is_primed());
    auto child_likelihoods = compute_likelihoods(child_genotypes, likelihood_model);
    if (debug_log_) debug::print(stream(*debug_log_), "child", child_likelihoods);
    boost::optional<double> lost_log_mass {};
    const auto reduced_child_likelihoods = reduce(child_likelihoods, prior_model_, lost_log_mass, options_);
    auto parent_likelihoods = compute_likelihoods(parent_genotypes, likelihood_model);
    if (debug_log_) debug::print(stream(*debug_log_), "parent", parent_likelihoods);
    const auto reduced_parent_likelihoods = reduce(parent_likelihoods, prior_model_, lost_log_mass, options_);
    haplotype_likelihoods.prime(trio_.child());
    auto joint_likelihoods = join(reduced_parent_likelihoods, reduced_child_likelihoods, mutation_model_);
    if (lost_log_mass) *lost_log_mass *= 2 * std::distance(reduced_child_likelihoods.first, reduced_child_likelihoods.last_to_join);
    clear(parent_likelihoods);
    clear(child_likelihoods);
    const auto evidence = normalise_exp(joint_likelihoods);
    if (debug_log_) debug::print(stream(*debug_log_), joint_likelihoods);
    return {std::move(joint_likelihoods), evidence, lost_log_mass};
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
