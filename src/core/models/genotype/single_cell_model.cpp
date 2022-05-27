// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "single_cell_model.hpp"

#include <utility>
#include <cassert>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <limits>

#include "utils/k_medoids.hpp"
#include "utils/select_top_k.hpp"
#include "utils/maths.hpp"
#include "utils/concat.hpp"
#include "subclone_model.hpp"
#include "population_model.hpp"
#include "coalescent_population_prior_model.hpp"
#include "genotype_prior_model.hpp"
#include "individual_model.hpp"

namespace octopus { namespace model {

const UniformPopulationPriorModel SingleCellModel::default_population_prior_model_{};

namespace {

auto get_vb_model_options(const SingleCellModel::AlgorithmParameters& config)
{
    VariationalBayesMixtureMixtureModel::Options result {};
    result.parallel_execution = config.execution_policy == ExecutionPolicy::par;
    return result;
}

} // namespace

SingleCellModel::SingleCellModel(std::vector<SampleName> samples,
                                 SingleCellPriorModel prior_model,
                                 Parameters parameters,
                                 AlgorithmParameters config,
                                 boost::optional<const PopulationPriorModel&> population_prior_model)
: samples_{std::move(samples)}
, prior_model_{std::move(prior_model)}
, posterior_model_ {get_vb_model_options(config)}
, parameters_{std::move(parameters)}
, config_{std::move(config)}
, population_prior_model_{std::addressof(population_prior_model ? *population_prior_model : default_population_prior_model_)}
{}

const SingleCellPriorModel& SingleCellModel::prior_model() const
{
    return prior_model_;
}

namespace {

template <typename Range>
bool all_same_ploidy(const Range& genotypes)
{
    const static auto ploidy_not_equal = [] (const auto& lhs, const auto& rhs) { return lhs.ploidy() != rhs.ploidy(); };
    return std::adjacent_find(std::cbegin(genotypes), std::cend(genotypes), ploidy_not_equal) == std::cend(genotypes);
}

template <typename T>
std::vector<T>
copy(const std::vector<T>& values, const std::vector<std::size_t>& indices)
{
    std::vector<T> result {};
    result.reserve(indices.size());
    for (auto idx : indices) {
        assert(idx < values.size());
        result.push_back(values[idx]);
    }
    return result;
}

} // namespace

SingleCellModel::Inferences
SingleCellModel::evaluate(const GenotypeVector& genotypes,
                          const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    assert(all_same_ploidy(genotypes));
    const auto num_clones = prior_model_.phylogeny().size();
    assert(num_clones <= genotypes.size());
    assert(num_clones <= samples_.size());
    Inferences result {};
    if (num_clones == 1) {
        evaluate(result, genotypes, haplotype_likelihoods);
    } else {
        const auto genotype_combinations = propose_genotype_combinations(genotypes, haplotype_likelihoods);
        evaluate(result, genotypes, genotype_combinations, haplotype_likelihoods);
    }
    return result;
}

SingleCellModel::Inferences
SingleCellModel::evaluate(const PhylogenyNodePloidyMap& phylogeny_ploidies,
                          const GenotypeVector& genotypes,
                          const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    assert(phylogeny_ploidies.size() == prior_model_.phylogeny().size());
    const auto num_clones = prior_model_.phylogeny().size();
    assert(num_clones <= genotypes.size());
    if (num_clones == 1) {
        return evaluate(genotypes, haplotype_likelihoods);
    } else {
        Inferences result {};
        const auto genotype_combinations = propose_genotype_combinations(phylogeny_ploidies, genotypes, haplotype_likelihoods);
        evaluate(result, genotypes, genotype_combinations, haplotype_likelihoods);
        return result;
    }
}

// private methods

std::vector<std::size_t>
SingleCellModel::propose_genotypes(const GenotypeVector& genotypes,
                                   const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    assert(config_.max_genotype_combinations);
    const auto merged_likelihoods = haplotype_likelihoods.merge_samples();
    const IndividualModel pooled_model {prior_model_.germline_prior_model()};
    const auto pooled_model_inferences = pooled_model.evaluate(genotypes, merged_likelihoods);
    return select_top_k_indices(pooled_model_inferences.posteriors.genotype_log_probabilities, *config_.max_genotype_combinations);
}

namespace {

auto log(std::size_t base, std::size_t x)
{
    return std::log2(x) / std::log2(base);
}

auto num_combinations(const std::size_t num_genotypes, const std::size_t num_samples)
{
    static constexpr auto max_combinations = std::numeric_limits<std::size_t>::max();
    if (num_samples <= log(num_genotypes, max_combinations)) {
        return static_cast<std::size_t>(std::pow(num_genotypes, num_samples));
    } else {
        return max_combinations;
    }
}

template <typename T>
auto select(const std::vector<std::size_t>& indices, const std::vector<T>& data)
{
    std::vector<T> result {};
    result.reserve(indices.size());
    for (auto idx : indices) result.push_back(data[idx]);
    return result;
}

auto l1_norm(const std::vector<double>& p, const std::vector<double>& q) noexcept
{
    return std::inner_product(std::cbegin(p), std::cend(p), std::cbegin(q), 0.0,
                              std::plus<> {}, [] (const auto a, const auto b) {
        return std::abs(a - b);
    });
}

ClusterVector
cluster_samples(const std::vector<PopulationModel::Latents::ProbabilityVector>& genotype_posteriors,
                const unsigned num_clusters)
{
    ClusterVector result {};
    KMediodsParameters params {};
    params.initialisation = KMediodsParameters::InitialisationMode::max_distance;
    auto best_fit = k_medoids(genotype_posteriors, num_clusters, result, l1_norm, params).second;
    params.initialisation = KMediodsParameters::InitialisationMode::total_distance;
    ClusterVector tmp {};
    auto fit = k_medoids(genotype_posteriors, num_clusters, tmp, l1_norm, params).second;
    if (fit < best_fit) {
        result = std::move(tmp);
        best_fit = fit;
    }
    params.initialisation = KMediodsParameters::InitialisationMode::random;
    for (int i {0}; i < 3; ++i) {
        tmp.clear();
        fit = k_medoids(genotype_posteriors, num_clusters, tmp, l1_norm, params).second;
        if (fit < best_fit) {
            result = std::move(tmp);
            best_fit = fit;
        }
    }
    return result;
}

template <typename ForwardIterator>
ForwardIterator
unique_stable(const ForwardIterator first, const ForwardIterator last)
{
    const auto n = static_cast<std::size_t>(std::distance(first, last));
    using ValueTp = typename std::iterator_traits<ForwardIterator>::value_type;
    std::vector<std::pair<ValueTp, std::size_t>> indexed_values(n);
    std::size_t idx {0};
    std::transform(std::make_move_iterator(first), std::make_move_iterator(last),
                   std::begin(indexed_values),
                   [&] (auto&& value) { return std::make_pair(std::move(value), idx++); });
    std::sort(std::begin(indexed_values), std::end(indexed_values));
    const static auto first_equal = [] (const auto& lhs, const auto& rhs) { return lhs.first == rhs.first; };
    indexed_values.erase(std::unique(std::begin(indexed_values), std::end(indexed_values), first_equal), std::end(indexed_values));
    const static auto index_less = [] (const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; };
    std::sort(std::begin(indexed_values), std::end(indexed_values), index_less);
    return std::transform(std::make_move_iterator(std::begin(indexed_values)),
                          std::make_move_iterator(std::end(indexed_values)),
                          first, [] (auto&& p) { return std::move(p.first); });
}

template <typename Range>
void unique_stable_erase(Range& values)
{
    values.erase(unique_stable(std::begin(values), std::end(values)), std::end(values));
}

template <typename T1, typename T2>
auto zip(std::vector<T1>&& lhs, std::vector<T2>&& rhs)
{
    assert(lhs.size() == rhs.size());
    std::vector<std::pair<T1, T2>> result {};
    result.reserve(lhs.size());
    std::transform(std::make_move_iterator(std::begin(lhs)), std::make_move_iterator(std::end(lhs)),
                   std::make_move_iterator(std::begin(rhs)), std::back_inserter(result),
                   [] (T1&& a, T2&& b) noexcept { return std::make_pair(std::move(a), std::move(b)); });
    return result;
}

template <typename T1, typename T2>
auto unzip(std::vector<std::pair<T1, T2>>&& zipped)
{
    std::vector<T1> lhs {}; std::vector<T2> rhs {};
    lhs.reserve(zipped.size()); rhs.reserve(zipped.size());
    for (auto& p : zipped) {
        lhs.push_back(std::move(p.first));
        rhs.push_back(std::move(p.second));
    }
    return std::make_pair(std::move(lhs), std::move(rhs));
}

void combine(const IndividualModel::Latents::ProbabilityVector& src, IndividualModel::Latents::ProbabilityVector& dst,
             const double src_weight = 1, const double dst_weight = 1) noexcept
{
    assert(src.size() == dst.size());
    const auto lhs_prob = src_weight / (src_weight + dst_weight);
    const auto rhs_prob = dst_weight / (src_weight + dst_weight);
    const auto combine_probs = [=] (auto lhs, auto rhs) noexcept { return lhs_prob * lhs + rhs_prob * rhs; };
    std::transform(std::cbegin(src), std::cend(src), std::cbegin(dst), std::begin(dst), combine_probs);
}

template <typename T>
std::vector<T> select(const std::vector<bool>& selectors, const std::vector<T>& values)
{
    assert(selectors.size() == values.size());
    std::vector<T> result {};
    result.reserve(std::accumulate(std::cbegin(selectors), std::cend(selectors), 0));
    std::size_t selector_idx {0};
    const auto selector = [&] (const auto&) { return selectors[selector_idx++]; };
    std::copy_if(std::cbegin(values), std::cend(values), std::back_inserter(result), selector);
    return result;
}

template <typename T>
IndexTupleVector
select_top_k_combinations(const std::vector<std::vector<T>>& values, const std::size_t k, const std::size_t n)
{
    const auto tuples = select_top_k_tuples(values, k);
    if (values.size() <= n) {
        return tuples;
    } else {
        IndexTupleVector result {};
        result.reserve(k);
        std::vector<bool> selectors(values.size());
        std::fill_n(std::rbegin(selectors), n, 1);
        for (const auto& tuple : tuples) {
            do {
                if (result.size() < k) {
                    result.push_back(select(selectors, tuple));
                } else {
                    return result;
                }
            } while (std::next_permutation(std::begin(selectors), std::end(selectors)));
        }
        return result;
    }
}

template <typename Range>
auto sum_entropies(const std::vector<Range>& ranges)
{
    const static auto add_entropy = [] (auto total, const auto& probs) { return maths::entropy2(probs); };
    return std::accumulate(std::cbegin(ranges), std::cend(ranges), 0.0, add_entropy);
}

template <typename Range>
auto max_entropy_idx(const std::vector<Range>& ranges)
{
    std::size_t result {0};
    double max {0};
    for (std::size_t idx {0}; idx < ranges.size(); ++idx) {
        const auto entropy = maths::entropy(ranges[idx]);
        if (idx == 0 || max < entropy) {
            result = idx;
            max = entropy;
        }
    }
    return result;
}

template <typename T>
void erase(std::vector<T>& values, const std::size_t idx)
{
    values.erase(std::next(std::cbegin(values), idx));
}

struct ClusterHelper
{
    std::vector<std::size_t> sample_indices;
    std::vector<SampleName> samples;
    ProbabilityVector genotype_probabilities;
    double entropy;
};

void sort_by_size_and_entropy(ClusterVector& clusters,
                              std::vector<std::vector<SampleName>>& samples_by_cluster,
                              std::vector<ProbabilityVector>& cluster_marginal_genotype_posteriors)
{
    std::vector<ClusterHelper> tmp_clusters(clusters.size());
    for (std::size_t k {0}; k < clusters.size(); ++k) {
        std::swap(tmp_clusters[k].sample_indices, clusters[k]);
        std::swap(tmp_clusters[k].samples, samples_by_cluster[k]);
        std::swap(tmp_clusters[k].genotype_probabilities, cluster_marginal_genotype_posteriors[k]);
        tmp_clusters[k].entropy = maths::entropy2(tmp_clusters[k].genotype_probabilities);
    }
    const static auto sorter = [] (const ClusterHelper& lhs, const ClusterHelper& rhs) {
        if (lhs.samples.size() == rhs.samples.size()) {
            return lhs.entropy < rhs.entropy;
        } else {
            return lhs.samples.size() > rhs.samples.size();
        }
    };
    std::sort(std::begin(tmp_clusters), std::end(tmp_clusters), sorter);
    for (std::size_t k {0}; k < clusters.size(); ++k) {
        std::swap(tmp_clusters[k].sample_indices, clusters[k]);
        std::swap(tmp_clusters[k].samples, samples_by_cluster[k]);
        std::swap(tmp_clusters[k].genotype_probabilities, cluster_marginal_genotype_posteriors[k]);
    }
}

void erase_combinations_with_duplicate_indices(std::vector<std::vector<std::size_t>>& combinations, const std::size_t max_index)
{
    std::vector<int> counts(max_index);
    combinations.erase(std::remove_if(std::begin(combinations), std::end(combinations), [&counts] (const auto& indices) {
        std::fill(std::begin(counts), std::end(counts), 0);
        for (auto idx : indices) ++counts[idx];
        return std::any_of(std::cbegin(counts), std::cend(counts), [] (auto count) { return count > 1; });
    }), std::end(combinations));
}

} // namespace

SingleCellModel::GenotypeCombinationVector
SingleCellModel::propose_genotype_combinations(const GenotypeVector& genotypes,
                                               const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    const auto num_groups = prior_model_.phylogeny().size();
    const auto max_possible_combinations = num_combinations(genotypes.size(), num_groups);
    const auto max_genotype_combinations = config_.max_genotype_combinations ? *config_.max_genotype_combinations : max_possible_combinations;
    
    // 1. Run population model
    // 2. Cluster samples
    // 3. Run individual model on merged reads
    // 4. Select top combinations using cluster marginal posteriors
    
    PopulationModel::Options population_model_options {};
    population_model_options.max_genotype_combinations = max_genotype_combinations;
    PopulationModel population_model {*population_prior_model_, population_model_options};
    std::vector<PopulationModel::Latents::ProbabilityVector> population_genotype_posteriors;
    IndividualModel individual_model {prior_model_.germline_prior_model()};
    std::vector<ProbabilityVector> cluster_marginal_genotype_posteriors {};
    cluster_marginal_genotype_posteriors.reserve(samples_.size() / 2);
    const auto haplotypes = haplotype_likelihoods.haplotypes();
    ClusterVector clusters {};
    std::vector<std::vector<SampleName>> samples_by_cluster {};
    
    while (clusters.empty() || clusters.size() > num_groups) {
        if (clusters.empty()) {
            auto population_inferences = population_model.evaluate(samples_, haplotypes, genotypes, haplotype_likelihoods);
            population_genotype_posteriors = std::move(population_inferences.posteriors.marginal_genotype_probabilities);
            clusters = cluster_samples(population_genotype_posteriors, std::min(std::max(samples_.size() / 4, 2 * num_groups), samples_.size()));
        } else if (clusters.size() > 2 * num_groups) {
            clusters = cluster_samples(cluster_marginal_genotype_posteriors, std::min(std::max(clusters.size() / 2, 2 * num_groups), samples_.size()));
        } else {
            auto num_effective_clusters = sum_entropies(cluster_marginal_genotype_posteriors);
            if ((clusters.size() == num_groups + 1 && num_effective_clusters < num_groups / 2)
             || (num_groups < 4 && clusters.size() <= num_groups + 2 && num_effective_clusters < 2)) {
                break;
            }
            sort_by_size_and_entropy(clusters, samples_by_cluster, cluster_marginal_genotype_posteriors);
            while (clusters.size() > num_groups && clusters.back().size() == 1
                && maths::entropy2(cluster_marginal_genotype_posteriors.back()) > 1) {
                clusters.pop_back();
                samples_by_cluster.pop_back();
                cluster_marginal_genotype_posteriors.pop_back();
            }
            if (clusters.size() <= num_groups + 1) break;
            num_effective_clusters = sum_entropies(cluster_marginal_genotype_posteriors);
            if (num_groups < 4 && clusters.size() <= num_groups + 2 && num_effective_clusters < 2) break;
            clusters = cluster_samples(cluster_marginal_genotype_posteriors, std::min(num_groups + 1, samples_.size()));
        }
        cluster_marginal_genotype_posteriors.clear();
        std::vector<std::vector<SampleName>> next_samples_by_cluster {};
        next_samples_by_cluster.reserve(clusters.size());
        for (const auto& cluster : clusters) {
            if (samples_by_cluster.empty()) {
                next_samples_by_cluster.push_back(select(cluster, samples_));
            } else {
                next_samples_by_cluster.push_back(concat(select(cluster, samples_by_cluster)));
            }
            const auto pooled_likelihoods = haplotype_likelihoods.merge_samples(next_samples_by_cluster.back());
            auto cluster_inferences = individual_model.evaluate(genotypes, pooled_likelihoods);
            cluster_marginal_genotype_posteriors.push_back(std::move(cluster_inferences.posteriors.genotype_probabilities));
        }
        samples_by_cluster = std::move(next_samples_by_cluster);
    }
    
    GenotypeCombinationVector result {};
    for (auto k = 10; result.empty() && k < 10'000; k *= 2) {
        result = select_top_k_combinations(cluster_marginal_genotype_posteriors, k * max_genotype_combinations, num_groups);
        // Remove combinations with duplicate genotypes as these are redundant according to model.
        erase_combinations_with_duplicate_indices(result, genotypes.size());
        // Reorder the combinations, keeping only the most probable one under the prior
        std::vector<SingleCellPriorModel::GenotypeReference> combination_refs {};
        combination_refs.reserve(num_groups);
        for (GenotypeCombination& combination : result) {
            std::sort(std::begin(combination), std::end(combination));
            GenotypeCombination best_combination;
            auto max_prior = std::numeric_limits<double>::lowest();
            do {
                for (auto idx : combination) combination_refs.emplace_back(genotypes[idx]);
                const auto prior = prior_model_.evaluate(combination_refs);
                combination_refs.clear();
                if (prior > max_prior) {
                    best_combination = combination;
                    max_prior = prior;
                }
            } while (std::next_permutation(std::begin(combination), std::end(combination)));
            combination = std::move(best_combination);
        }
        unique_stable_erase(result);
    }
    if (result.size() > max_genotype_combinations) {
        result.resize(max_genotype_combinations);
    } else if (result.empty()) {
        GenotypeCombination combo(num_groups);
        std::iota(std::begin(combo), std::end(combo), 0);
        result.push_back(std::move(combo));
    }
    return result;
}

SingleCellModel::GenotypeCombinationVector
SingleCellModel::propose_all_genotype_combinations(const GenotypeVector& genotypes) const
{
    const auto num_groups = prior_model_.phylogeny().size();
    GenotypeCombinationVector result {};
    result.reserve(num_combinations(genotypes.size(), num_groups));
    GenotypeCombination tmp(num_groups);
    std::vector<bool> v(genotypes.size() * num_groups), hits(genotypes.size(), false);
    std::fill(std::begin(v), std::next(std::begin(v), num_groups), true);
    do {
        bool good {true};
        for (std::size_t i {0}, k {0}; k < num_groups; ++i) {
            if (v[i]) {
                if (i / genotypes.size() == k) {
                    tmp[k++] = i - genotypes.size() * (i / genotypes.size());
                } else {
                    good = false;
                    k = num_groups;
                }
            }
        }
        if (good) {
            for (auto idx : tmp) {
                if (hits[idx]) {
                    good = false;
                    break;
                } else {
                    hits[idx] = true;
                }
            }
            std::fill(std::begin(hits), std::end(hits), false);
            if (good) result.push_back(tmp);
        }
    } while (std::prev_permutation(std::begin(v), std::end(v)));
    return result;
}

namespace {

auto max_ploidy(const SingleCellModel::GenotypeVector& genotypes)
{
    const static auto ploidy_less = [] (const auto& lhs, const auto& rhs) { return lhs.ploidy() < rhs.ploidy(); };
    return std::max_element(std::cbegin(genotypes), std::cend(genotypes), ploidy_less)->ploidy();
}

auto segment_by_ploidy(const SingleCellModel::GenotypeVector& genotypes)
{
    std::vector<SingleCellModel::GenotypeVector> result(max_ploidy(genotypes) + 1);
    for (auto& g : result) g.reserve(genotypes.size());
    for (const auto& genotype : genotypes) result[genotype.ploidy()].push_back(genotype);
    return result;
}

bool valid_ploidies(const std::vector<std::size_t>& combination,
                    const SingleCellModel::GenotypeVector& genotypes,
                    const std::vector<unsigned>& required_ploidies) noexcept
{
    assert(combination.size() == required_ploidies.size());
    std::size_t idx {0};
    const auto ploidy_equal = [&] (auto g) noexcept { return genotypes[g].ploidy() == required_ploidies[idx++]; };
    return std::all_of(std::cbegin(combination), std::cend(combination), ploidy_equal);
}

template <typename T>
auto get_leaf_labels(const Phylogeny<std::size_t, T>& phylogeny)
{
    std::vector<std::size_t> result {};
    result.reserve(phylogeny.size() - 1);
    for (std::size_t id {0}; id < phylogeny.size(); ++id) {
        if (phylogeny.num_descendants(id) == 0) {
            result.push_back(id);
        }
    }
    return result;
}

bool is_redundant_copy_change(const Genotype<IndexedHaplotype<>>& lhs, const Genotype<IndexedHaplotype<>>& rhs)
{
    return lhs.ploidy() != rhs.ploidy() && have_same_elements(lhs, rhs);
}

void erase_combinations_with_redundant_copy_changes(std::vector<std::vector<std::size_t>>& combinations,
                                                    const SingleCellModel::GenotypeVector& genotypes,
                                                    const SingleCellPriorModel::CellPhylogeny& phylogeny)
{
    const auto leafs = get_leaf_labels(phylogeny);
    std::vector<std::size_t> leaf_ancestors(leafs.size());
    std::transform(std::cbegin(leafs), std::cend(leafs), std::begin(leaf_ancestors),
                   [&] (auto leaf_id) { return phylogeny.ancestor(leaf_id).id; });
    std::vector<std::vector<boost::optional<bool>>> have_same_elements_cache(genotypes.size(), std::vector<boost::optional<bool>>(genotypes.size()));
    combinations.erase(std::remove_if(std::begin(combinations), std::end(combinations), [&] (const auto& combination) {
        std::size_t leaf_idx {0};
        return std::any_of(std::cbegin(leafs), std::cend(leafs), [&] (const auto leaf) {
            const auto child_index = combination[leaf];
            const auto parent_index = combination[leaf_ancestors[leaf_idx++]];
            auto& result = have_same_elements_cache[child_index][parent_index];
            if (!result) result = is_redundant_copy_change(genotypes[child_index], genotypes[parent_index]);
            return *result;
        });
    }), std::end(combinations));
}

} // namespace

class ZygosityGenotypePriorModel : public GenotypePriorModel
{
public:
    using GenotypePriorModel::LogProbability;
    
    ZygosityGenotypePriorModel(const GenotypePriorModel& base) : base_ {std::addressof(base)} {}
    
    virtual ~ZygosityGenotypePriorModel() = default;

private:
    const GenotypePriorModel* base_;
    
    template <typename T>
    LogProbability do_evaluate_helper(const Genotype<T>& genotype) const
    {
        LogProbability result {(ploidy(genotype) - zygosity(genotype)) * std::log(0.1)};
        if (base_) {
            result += base_->evaluate(genotype);
        }
        return result;
    }
    virtual LogProbability do_evaluate(const Genotype<Haplotype>& genotype) const override { return do_evaluate_helper(genotype);}
    virtual LogProbability do_evaluate(const Genotype<IndexedHaplotype<>>& genotype) const override { return do_evaluate_helper(genotype); }
    bool check_is_primed() const noexcept override { return true; }
    
};

SingleCellModel::GenotypeCombinationVector
SingleCellModel::propose_genotype_combinations(const PhylogenyNodePloidyMap& phylogeny_ploidies,
                                               const GenotypeVector& genotypes,
                                               const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    // First assign each sample to a ploidy
    const auto genotypes_by_ploidy = segment_by_ploidy(genotypes);
    std::vector<unsigned> sample_ploidies {};
    sample_ploidies.reserve(samples_.size());
    const ZygosityGenotypePriorModel zygosity_prior {prior_model_.germline_prior_model()};
    const IndividualModel zygosity_individual_model {zygosity_prior};
    for (const auto& sample : samples_) {
        haplotype_likelihoods.prime(sample);
        unsigned best_ploidy {1};
        double max_log_evidence {};
        for (unsigned ploidy {1}; ploidy < genotypes_by_ploidy.size(); ++ploidy) {
            if (!genotypes_by_ploidy[ploidy].empty()) {
                const auto inferences = zygosity_individual_model.evaluate(genotypes_by_ploidy[ploidy], haplotype_likelihoods);
                if (ploidy == 1 || max_log_evidence < inferences.log_evidence) {
                    best_ploidy = ploidy;
                    max_log_evidence = inferences.log_evidence;
                }
            }
        }
        sample_ploidies.push_back(best_ploidy);
    }
    
    PopulationModel::Options population_model_options {};
    population_model_options.max_genotype_combinations = config_.max_genotype_combinations;
    PopulationModel population_model {*population_prior_model_, population_model_options};
    std::vector<PopulationModel::Latents::ProbabilityVector> population_genotype_posteriors;
    const auto num_groups = prior_model_.phylogeny().size();
    std::vector<ProbabilityVector> cluster_marginal_genotype_posteriors {};
    const auto haplotypes = haplotype_likelihoods.haplotypes();
    ClusterVector clusters {};
    std::vector<std::vector<SampleName>> samples_by_cluster {};
    
    while (clusters.empty() || clusters.size() > num_groups) {
        if (clusters.empty()) {
            auto population_inferences = population_model.evaluate(samples_, haplotypes, genotypes, haplotype_likelihoods);
            population_genotype_posteriors = std::move(population_inferences.posteriors.marginal_genotype_probabilities);
            clusters = cluster_samples(population_genotype_posteriors, std::min(std::max(samples_.size() / 4, 2 * num_groups), samples_.size()));
        } else if (clusters.size() > 2 * num_groups) {
            clusters = cluster_samples(cluster_marginal_genotype_posteriors, std::min(std::max(clusters.size() / 2, 2 * num_groups), samples_.size()));
        } else {
            auto num_effective_clusters = sum_entropies(cluster_marginal_genotype_posteriors);
            if ((clusters.size() == num_groups + 1 && num_effective_clusters < num_groups / 2)
                || (num_groups < 4 && clusters.size() <= num_groups + 2 && num_effective_clusters < 1)) {
                break;
            }
            sort_by_size_and_entropy(clusters, samples_by_cluster, cluster_marginal_genotype_posteriors);
            while (clusters.size() > num_groups && clusters.back().size() == 1
                   && maths::entropy2(cluster_marginal_genotype_posteriors.back()) > 1) {
                clusters.pop_back();
                samples_by_cluster.pop_back();
                cluster_marginal_genotype_posteriors.pop_back();
            }
            if (clusters.size() <= num_groups + 1) break;
            clusters = cluster_samples(cluster_marginal_genotype_posteriors, std::min(num_groups + 1, samples_.size()));
        }
        cluster_marginal_genotype_posteriors.clear();
        std::vector<std::vector<SampleName>> next_samples_by_cluster {};
        next_samples_by_cluster.reserve(clusters.size());
        for (const auto& cluster : clusters) {
            if (samples_by_cluster.empty()) {
                next_samples_by_cluster.push_back(select(cluster, samples_));
            } else {
                next_samples_by_cluster.push_back(concat(select(cluster, samples_by_cluster)));
            }
            const auto pooled_likelihoods = haplotype_likelihoods.merge_samples(next_samples_by_cluster.back());
            auto cluster_inferences = zygosity_individual_model.evaluate(genotypes, pooled_likelihoods);
            cluster_marginal_genotype_posteriors.push_back(std::move(cluster_inferences.posteriors.genotype_probabilities));
        }
        samples_by_cluster = std::move(next_samples_by_cluster);
    }
    
    std::vector<unsigned> required_ploidies(num_groups);
    for (std::size_t id {0}; id < num_groups; ++id) required_ploidies[id] = phylogeny_ploidies.at(id);
    GenotypeCombinationVector result {};
    const auto k = config_.max_genotype_combinations ? 100 * *config_.max_genotype_combinations : std::numeric_limits<std::size_t>::max();
    result = select_top_k_combinations(cluster_marginal_genotype_posteriors, k, num_groups);
    // Remove combinations with duplicate genotypes as these are redundant according to model.
    erase_combinations_with_duplicate_indices(result, genotypes.size());
    // Reorder the combinations, keeping only the most probable one under the prior
    std::vector<SingleCellPriorModel::GenotypeReference> combination_refs {};
    combination_refs.reserve(num_groups);
    result.erase(std::remove_if(std::begin(result), std::end(result), [&] (GenotypeCombination& combination) {
        std::sort(std::begin(combination), std::end(combination));
        GenotypeCombination best_combination;
        auto max_prior = std::numeric_limits<double>::lowest();
        bool has_valid_combination {false};
        do {
            if (valid_ploidies(combination, genotypes, required_ploidies)) {
                for (auto idx : combination) combination_refs.emplace_back(genotypes[idx]);
                const auto prior = prior_model_.evaluate(combination_refs);
                combination_refs.clear();
                if (prior > max_prior) {
                    best_combination = combination;
                    max_prior = prior;
                }
                has_valid_combination = true;
            }
        } while (std::next_permutation(std::begin(combination), std::end(combination)));
        combination = std::move(best_combination);
        return !has_valid_combination;
    }), std::end(result));
    unique_stable_erase(result);
    erase_combinations_with_redundant_copy_changes(result, genotypes, prior_model_.phylogeny());
    if (result.empty()) {
        throw NoViableGenotypeCombinationsError {};
    }
    if (config_.max_genotype_combinations && result.size() > *config_.max_genotype_combinations) {
        result.resize(*config_.max_genotype_combinations);
    }
    return result;
}

namespace {

auto make_point_seed(const std::size_t num_genotypes, const std::size_t idx, const double p = 0.9999)
{
    LogProbabilityVector result(num_genotypes, num_genotypes > 1 ? std::log((1 - p) / (num_genotypes - 1)) : 0);
    if (num_genotypes > 1) result[idx] = std::log(p);
    return result;
}

template <typename InputIterator>
void make_point_seeds(const std::size_t num_genotypes,
                      InputIterator first_indices_itr, InputIterator last_indices_itr,
                      std::vector<LogProbabilityVector>& result, const double p = 0.9999)
{
    result.reserve(result.size() + std::distance(first_indices_itr, last_indices_itr));
    std::transform(first_indices_itr, last_indices_itr, std::back_inserter(result),
                   [=] (auto idx) { return make_point_seed(num_genotypes, idx, p); });
}

void make_point_seeds(const std::size_t num_genotypes, const std::vector<std::size_t>& indices,
                      std::vector<LogProbabilityVector>& result, const double p = 0.9999)
{
    make_point_seeds(num_genotypes, std::cbegin(indices), std::cend(indices), result, p);
}

} // namespace

void
SingleCellModel::evaluate(Inferences& result,
                          const GenotypeVector& genotypes,
                          const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    SubcloneModel::Priors subclone_priors {prior_model_.germline_prior_model(), {}};
    const auto ploidy = genotypes.front().ploidy();
    for (const auto& sample : samples_) {
        subclone_priors.alphas.emplace(sample, SubcloneModel::Priors::GenotypeMixturesDirichletAlphas(ploidy,parameters_.dropout_concentration));
    }
    SubcloneModel::AlgorithmParameters algo_params {};
    algo_params.max_seeds = config_.max_seeds;
    SubcloneModel helper_model {samples_, std::move(subclone_priors), algo_params};
    SubcloneModel::InferredLatents subclone_inferences;
    if (!config_.max_genotype_combinations || genotypes.size() <= *config_.max_genotype_combinations) {
        subclone_inferences = helper_model.evaluate(genotypes, haplotype_likelihoods);
    } else {
        const auto genotype_subset_indices = propose_genotypes(genotypes, haplotype_likelihoods);
        const auto genotype_subset = copy(genotypes, genotype_subset_indices);
        const auto num_hints = std::min(genotype_subset_indices.size(), static_cast<std::size_t>(config_.max_seeds / 2));
        std::vector<std::size_t> hint_indices(num_hints);
        std::iota(std::begin(hint_indices), std::end(hint_indices), 0u);
        std::vector<LogProbabilityVector> hints {};
        make_point_seeds(genotype_subset.size(), hint_indices, hints);
        subclone_inferences = helper_model.evaluate(genotype_subset, haplotype_likelihoods, std::move(hints));
        SubcloneModel::Latents::ProbabilityVector weighted_genotype_posteriors(genotypes.size());
        for (std::size_t g {0}; g < genotype_subset.size(); ++g) {
            weighted_genotype_posteriors[genotype_subset_indices[g]] = subclone_inferences.weighted_genotype_posteriors[g];
        }
        subclone_inferences.weighted_genotype_posteriors = std::move(weighted_genotype_posteriors);
    }
    Inferences::GroupInferences founder {};
    founder.genotype_posteriors = std::move(subclone_inferences.weighted_genotype_posteriors);
    founder.sample_attachment_posteriors.assign(samples_.size(), 1.0);
    result.phylogeny.set_founder({0, std::move(founder)});
    result.log_evidence = subclone_inferences.approx_log_evidence;
}

void
SingleCellModel::evaluate(Inferences& result,
                          const GenotypeVector& genotypes,
                          const GenotypeCombinationVector& genotype_combinations,
                          const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    const auto genotype_combination_priors = calculate_genotype_priors(genotype_combinations, genotypes);
    const auto vb_haplotype_likelihoods = make_likelihood_matrix(genotype_combinations, genotypes, haplotype_likelihoods);
    auto seeds = propose_seeds(genotype_combinations, genotypes, genotype_combination_priors, haplotype_likelihoods);
    auto vb_inferences = evaluate_model(genotype_combination_priors, vb_haplotype_likelihoods, std::move(seeds));
    for (std::size_t group_idx {0}; group_idx < prior_model_.phylogeny().size(); ++group_idx) {
        Inferences::GroupInferences group {};
        group.sample_attachment_posteriors.resize(samples_.size());
        for (std::size_t sample_idx {0}; sample_idx < samples_.size(); ++sample_idx) {
            group.sample_attachment_posteriors[sample_idx] = vb_inferences.weighted_group_responsibilities[sample_idx][group_idx];
        }
        // Marginalise over genotypes
        group.genotype_posteriors.resize(genotypes.size());
        for (std::size_t genotype_combo_idx {0};
             genotype_combo_idx < genotype_combinations.size(); ++genotype_combo_idx) {
            group.genotype_posteriors[genotype_combinations[genotype_combo_idx][group_idx]] += vb_inferences.weighted_genotype_posteriors[genotype_combo_idx];
        }
        if (group_idx == 0) {
            result.phylogeny.set_founder({group_idx, std::move(group)});
        } else {
            const auto ancestor_idx = prior_model_.phylogeny().ancestor(group_idx).id;
            result.phylogeny.add_descendant({group_idx, std::move(group)}, ancestor_idx);
        }
    }
    result.log_evidence = vb_inferences.map.approx_log_evidence;
}

VariationalBayesMixtureMixtureModel::LogProbabilityVector
SingleCellModel::calculate_genotype_priors(const GenotypeCombinationVector& genotype_combinations,
                                           const GenotypeVector& genotypes) const
{
    std::vector<SingleCellPriorModel::GenotypeReference> combination_buffer {};
    LogProbabilityVector result(genotype_combinations.size());
    std::transform(std::cbegin(genotype_combinations), std::cend(genotype_combinations), std::begin(result),
                   [&] (const auto& combination) {
                       combination_buffer.clear();
                       combination_buffer.reserve(combination.size());
                       std::transform(std::cbegin(combination), std::cend(combination), std::back_inserter(combination_buffer),
                                      [&] (auto idx) { return std::cref(genotypes[idx]); });
                       return prior_model_.evaluate(combination_buffer);
    });
    return result;
}

SingleCellModel::VBLikelihoodMatrix
SingleCellModel::make_likelihood_matrix(const GenotypeCombinationVector& genotype_combinations,
                                        const GenotypeVector& genotypes,
                                        const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    VBLikelihoodMatrix result {};
    result.reserve(samples_.size());
    for (const auto& sample : samples_) {
        haplotype_likelihoods.prime(sample);
        VariationalBayesMixtureMixtureModel::GenotypeCombinationLikelihoodVector vb_combination_likelihoods {};
        vb_combination_likelihoods.reserve(genotype_combinations.size());
        for (const auto& genotype_combination : genotype_combinations) {
            VariationalBayesMixtureMixtureModel::GenotypeLikelihoodVector vb_genotype_likelihoods {};
            vb_genotype_likelihoods.reserve(genotype_combination.size());
            for (const auto genotype_idx : genotype_combination) {
                VariationalBayesMixtureMixtureModel::HaplotypeLikelihoodVector vb_haplotype_likelihoods {};
                vb_haplotype_likelihoods.reserve(genotypes[genotype_idx].ploidy());
                for (const auto& haplotype : genotypes[genotype_idx]) {
                    vb_haplotype_likelihoods.emplace_back(haplotype_likelihoods[haplotype]);
                }
                vb_genotype_likelihoods.push_back(std::move(vb_haplotype_likelihoods));
            }
            vb_combination_likelihoods.push_back(std::move(vb_genotype_likelihoods));
        }
        result.push_back(std::move(vb_combination_likelihoods));
    }
    return result;
}

namespace {

LogProbabilityVector log_uniform_dist(const std::size_t n)
{
    return LogProbabilityVector(n, -std::log(static_cast<double>(n)));
}

auto make_range_seed(const std::size_t num_genotypes, const std::size_t begin, const std::size_t n,
                     const double p = 0.9999999)
{
    LogProbabilityVector result(num_genotypes, std::log((1 - p) / (num_genotypes - n)));
    std::fill_n(std::next(std::begin(result), begin), n, std::log(p / n));
    return result;
}

} // namespace

SingleCellModel::VBSeedVector
SingleCellModel::propose_seeds(const GenotypeCombinationVector& genotype_combinations,
                               const GenotypeVector& genotypes,
                               const VariationalBayesMixtureMixtureModel::LogProbabilityVector& genotype_combination_priors,
                               const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    VBSeedVector result {};
    const auto num_genotypes = genotype_combinations.size();
    result.push_back(log_uniform_dist(num_genotypes));
    if (config_.max_seeds > 1 && num_genotypes > 100) result.push_back(make_range_seed(num_genotypes, 0, num_genotypes / 100));
    if (config_.max_seeds > 2 && num_genotypes > 10) result.push_back(make_range_seed(num_genotypes, 0, num_genotypes / 10));
    std::vector<std::size_t> top_indices(std::min(config_.max_seeds - result.size(), num_genotypes));
    std::iota(std::begin(top_indices), std::end(top_indices), 0u);
    make_point_seeds(genotype_combinations.size(), top_indices, result);
    return result;
}

VariationalBayesMixtureMixtureModel::Inferences
SingleCellModel::evaluate_model(const VariationalBayesMixtureMixtureModel::LogProbabilityVector& genotype_combination_priors,
                                const VBLikelihoodMatrix& haplotype_likelihoods,
                                VBSeedVector seeds) const
{
    if (parameters_.group_priors) {
        if (parameters_.sample_dropout_concentrations.size() == samples_.size()) {
            return posterior_model_.evaluate(genotype_combination_priors,
                                             haplotype_likelihoods,
                                             *parameters_.group_priors,
                                             parameters_.group_concentration,
                                             parameters_.sample_dropout_concentrations,
                                             std::move(seeds));
        } else {
            return posterior_model_.evaluate(genotype_combination_priors,
                                             haplotype_likelihoods,
                                             *parameters_.group_priors,
                                             parameters_.group_concentration,
                                             parameters_.dropout_concentration,
                                             std::move(seeds));
        }
    } else {
        if (parameters_.sample_dropout_concentrations.size() == samples_.size()) {
            return posterior_model_.evaluate(genotype_combination_priors,
                                             haplotype_likelihoods,
                                             parameters_.group_concentration,
                                             parameters_.sample_dropout_concentrations,
                                             std::move(seeds));
        } else {
            return posterior_model_.evaluate(genotype_combination_priors,
                                             haplotype_likelihoods,
                                             parameters_.group_concentration,
                                             parameters_.dropout_concentration,
                                             std::move(seeds));
        }
    }
}

} // namespace model
} // namespace octopus
