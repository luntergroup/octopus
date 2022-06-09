// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef subclone_model_hpp
#define subclone_model_hpp

#include <vector>
#include <unordered_map>
#include <utility>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/indexed_haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/types/cancer_genotype.hpp"
#include "core/types/partitioned_genotype.hpp"
#include "core/models/haplotype_likelihood_array.hpp"
#include "containers/mappable_block.hpp"
#include "exceptions/unimplemented_feature_error.hpp"
#include "utils/thread_pool.hpp"
#include "variational_bayes_mixture_model.hpp"
#include "genotype_prior_model.hpp"
#include "cancer_genotype_prior_model.hpp"
#include "haplogroup_prior_model.hpp"

namespace octopus { namespace model {

template <typename GenotypeType, typename GenotypePriorModel_>
class SubcloneModelBase
{
public:
    constexpr static unsigned max_ploidy = 100;
    struct AlgorithmParameters
    {
        unsigned max_iterations = 1000;
        double epsilon          = 0.05;
        unsigned max_seeds      = 12;
        boost::optional<MemoryFootprint> target_max_memory = boost::none;
    };
    
    struct Priors
    {
        using GenotypeMixturesDirichletAlphas   = std::vector<double>;
        using GenotypeMixturesDirichletAlphaMap = std::unordered_map<SampleName, GenotypeMixturesDirichletAlphas>;
        
        const GenotypePriorModel_& genotype_prior_model;
        GenotypeMixturesDirichletAlphaMap alphas;
    };
    
    struct Latents
    {
        using GenotypeMixturesDirichletAlphas   = std::vector<double>;
        using GenotypeMixturesDirichletAlphaMap = std::unordered_map<SampleName, GenotypeMixturesDirichletAlphas>;
        using ProbabilityVector                 = std::vector<double>;
        using LogProbabilityVector              = std::vector<double>;
        
        ProbabilityVector genotype_probabilities;
        LogProbabilityVector genotype_log_probabilities;
        GenotypeMixturesDirichletAlphaMap alphas;
    };
    
    struct InferredLatents
    {
        Latents max_evidence_params;
        typename Latents::LogProbabilityVector genotype_log_priors;
        typename Latents::ProbabilityVector weighted_genotype_posteriors;
        double approx_log_evidence;
    };

    using OptionalThreadPool = boost::optional<ThreadPool&>;
    
    SubcloneModelBase() = delete;
    
    SubcloneModelBase(std::vector<SampleName> samples, Priors priors);
    SubcloneModelBase(std::vector<SampleName> samples, Priors priors, AlgorithmParameters parameters);
    
    SubcloneModelBase(const SubcloneModelBase&)            = default;
    SubcloneModelBase& operator=(const SubcloneModelBase&) = default;
    SubcloneModelBase(SubcloneModelBase&&)                 = default;
    SubcloneModelBase& operator=(SubcloneModelBase&&)      = default;
    
    ~SubcloneModelBase() = default;
    
    const Priors& priors() const noexcept;
    
    void prime(const MappableBlock<Haplotype>& haplotypes);
    void unprime() noexcept;
    bool is_primed() const noexcept;
    
    InferredLatents
    evaluate(const std::vector<GenotypeType>& genotypes,
             const HaplotypeLikelihoodArray& haplotype_likelihoods,
             std::vector<typename Latents::LogProbabilityVector> hints = {},
             OptionalThreadPool = boost::none) const;
    
private:
    std::vector<SampleName> samples_;
    Priors priors_;
    AlgorithmParameters parameters_;
    const MappableBlock<Haplotype>* haplotypes_;
};

using SubcloneModel = SubcloneModelBase<Genotype<IndexedHaplotype<>>, GenotypePriorModel>;
using SomaticSubcloneModel = SubcloneModelBase<CancerGenotype<IndexedHaplotype<>>, CancerGenotypePriorModel>;
using HaplogroupSubcloneModel = SubcloneModelBase<PartitionedGenotype<IndexedHaplotype<>>, HaplogroupGenotypePriorModel>;

template <typename G, typename GPM>
constexpr unsigned SubcloneModelBase<G, GPM>::max_ploidy;

template <typename G, typename GPM>
SubcloneModelBase<G, GPM>::SubcloneModelBase(std::vector<SampleName> samples, Priors priors)
: SubcloneModelBase {std::move(samples), std::move(priors), AlgorithmParameters {}}
{}

template <typename G, typename GPM>
SubcloneModelBase<G, GPM>::SubcloneModelBase(std::vector<SampleName> samples, Priors priors, AlgorithmParameters parameters)
: samples_ {std::move(samples)}
, priors_ {std::move(priors)}
, parameters_ {parameters}
{}

template <typename G, typename GPM>
const typename SubcloneModelBase<G, GPM>::Priors& SubcloneModelBase<G, GPM>::priors() const noexcept
{
    return priors_;
}

template <typename G, typename GPM>
void SubcloneModelBase<G, GPM>::prime(const MappableBlock<Haplotype>& haplotypes)
{
    haplotypes_ = std::addressof(haplotypes);
}

template <typename G, typename GPM>
void SubcloneModelBase<G, GPM>::unprime() noexcept
{
    haplotypes_ = nullptr;
}

template <typename G, typename GPM>
bool SubcloneModelBase<G, GPM>::is_primed() const noexcept
{
    return haplotypes_;
}

namespace detail {

std::vector<LogProbabilityVector>
generate_seeds(const std::vector<SampleName>& samples,
               const std::vector<Genotype<IndexedHaplotype<>>>& genotypes,
               const LogProbabilityVector& genotype_log_priors,
               const HaplotypeLikelihoodArray& haplotype_log_likelihoods,
               const SubcloneModel::Priors& priors,
               std::size_t max_seeds,
               std::vector<LogProbabilityVector> hints = {},
               const MappableBlock<Haplotype>* haplotypes = nullptr);

std::vector<LogProbabilityVector>
generate_seeds(const std::vector<SampleName>& samples,
               const std::vector<CancerGenotype<IndexedHaplotype<>>>& genotypes,
               const LogProbabilityVector& genotype_log_priors,
               const HaplotypeLikelihoodArray& haplotype_log_likelihoods,
               const SomaticSubcloneModel::Priors& priors,
               std::size_t max_seeds,
               std::vector<LogProbabilityVector> hints = {},
               const MappableBlock<Haplotype>* haplotypes = nullptr);

std::vector<LogProbabilityVector>
generate_seeds(const std::vector<SampleName>& samples,
               const std::vector<PartitionedGenotype<IndexedHaplotype<>>>& genotypes,
               const LogProbabilityVector& genotype_log_priors,
               const HaplotypeLikelihoodArray& haplotype_log_likelihoods,
               const HaplogroupSubcloneModel::Priors& priors,
               std::size_t max_seeds,
               std::vector<LogProbabilityVector> hints = {},
               const MappableBlock<Haplotype>* haplotypes = nullptr);

template <std::size_t K, typename G, typename GPM>
VBAlpha<K> flatten(const typename SubcloneModelBase<G, GPM>::Priors::GenotypeMixturesDirichletAlphas& alpha)
{
    VBAlpha<K> result {};
    std::copy_n(std::cbegin(alpha), K, std::begin(result));
    return result;
}

template <std::size_t K, typename G, typename GPM>
VBAlphaVector<K> flatten(const typename SubcloneModelBase<G, GPM>::Priors::GenotypeMixturesDirichletAlphaMap& alphas,
                         const std::vector<SampleName>& samples)
{
    VBAlphaVector<K> result(samples.size());
    std::transform(std::cbegin(samples), std::cend(samples), std::begin(result),
                   [&alphas] (const auto& sample) { return flatten<K, G, GPM>(alphas.at(sample)); });
    return result;
}

template <std::size_t K>
VBGenotype<K>
flatten(const Genotype<IndexedHaplotype<>>& genotype,
        const SampleName& sample,
        const HaplotypeLikelihoodArray& haplotype_likelihoods)
{
    VBGenotype<K> result {};
    haplotype_likelihoods.prime(sample);
    std::transform(std::cbegin(genotype), std::cend(genotype), std::begin(result),
                   [&] (const auto& haplotype) { return std::cref(haplotype_likelihoods[haplotype]); });
    return result;
}

template <std::size_t K>
auto copy_cref(const Genotype<IndexedHaplotype<>>& genotype,
               const SampleName& sample,
               const HaplotypeLikelihoodArray& haplotype_likelihoods,
               typename VBGenotype<K>::iterator result_itr)
{
    haplotype_likelihoods.prime(sample);
    return std::transform(std::cbegin(genotype), std::cend(genotype), result_itr,
                          [&] (const auto& haplotype) { return std::cref(haplotype_likelihoods[haplotype]); });
}

template <std::size_t K>
VBGenotype<K>
flatten(const CancerGenotype<IndexedHaplotype<>>& genotype,
        const SampleName& sample,
        const HaplotypeLikelihoodArray& haplotype_likelihoods)
{
    VBGenotype<K> result {};
    assert(genotype.ploidy() == K);
    auto itr = copy_cref<K>(genotype.germline(), sample, haplotype_likelihoods, std::begin(result));
    copy_cref<K>(genotype.somatic(), sample, haplotype_likelihoods, itr);
    return result;
}

template <std::size_t K>
VBGenotype<K>
flatten(const PartitionedGenotype<IndexedHaplotype<>>& genotype,
        const SampleName& sample,
        const HaplotypeLikelihoodArray& haplotype_likelihoods)
{
    VBGenotype<K> result {};
    assert(genotype.ploidy() == K);
    auto itr = std::begin(result);
    for (std::size_t n {0}; n < genotype.num_partitions(); ++n) {
        itr = copy_cref<K>(genotype.partition(n), sample, haplotype_likelihoods, itr);
    }
    return result;
}

template <std::size_t K, typename G>
VBGenotypeVector<K>
flatten(const std::vector<G>& genotypes, const SampleName& sample,
        const HaplotypeLikelihoodArray& haplotype_likelihoods)
{
    VBGenotypeVector<K> result(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                   [&sample, &haplotype_likelihoods] (const auto& genotype) {
                       return flatten<K>(genotype, sample, haplotype_likelihoods);
                   });
    return result;
}

template <std::size_t K, typename G>
VBReadLikelihoodMatrix<K>
flatten(const std::vector<G>& genotypes,
        const std::vector<SampleName>& samples,
        const HaplotypeLikelihoodArray& haplotype_likelihoods)
{
    VBReadLikelihoodMatrix<K> result {};
    result.reserve(samples.size());
    std::transform(std::cbegin(samples), std::cend(samples), std::back_inserter(result),
                   [&genotypes, &haplotype_likelihoods] (const auto& sample) {
                       return flatten<K>(genotypes, sample, haplotype_likelihoods);
                   });
    return result;
}

template <std::size_t K, typename G, typename GPM>
auto expand(VBAlpha<K>& alpha)
{
    return typename SubcloneModelBase<G, GPM>::Latents::GenotypeMixturesDirichletAlphas(std::begin(alpha), std::end(alpha));
}

template <std::size_t K, typename G, typename GPM>
auto expand(const std::vector<SampleName>& samples, VBAlphaVector<K>&& alphas)
{
    typename SubcloneModelBase<G, GPM>::Latents::GenotypeMixturesDirichletAlphaMap result {};
    std::transform(std::cbegin(samples), std::cend(samples), std::begin(alphas),
                   std::inserter(result, std::begin(result)),
                   [] (const auto& sample, auto&& vb_alpha) {
                       return std::make_pair(sample, expand<K, G, GPM>(vb_alpha));
                   });
    return result;
}

template <std::size_t K, typename G, typename GPM>
typename SubcloneModelBase<G, GPM>::InferredLatents
expand(VBResultPacket<K>& vb_results, const std::vector<SampleName>& samples, LogProbabilityVector genotype_log_priors)
{
    typename SubcloneModelBase<G, GPM>::Latents posterior_latents {
        std::move(vb_results.map_latents.genotype_posteriors),
        std::move(vb_results.map_latents.genotype_log_posteriors),
        expand<K, G, GPM>(samples, std::move(vb_results.map_latents.alphas))
        };
    return {std::move(posterior_latents),
            std::move(genotype_log_priors),
            std::move(vb_results.evidence_weighted_genotype_posteriors),
            vb_results.max_log_evidence};
}

template <std::size_t K, typename G, typename GPM>
typename SubcloneModelBase<G, GPM>::InferredLatents
run_variational_bayes_helper(const std::vector<SampleName>& samples,
                             const std::vector<G>& genotypes,
                             const typename SubcloneModelBase<G, GPM>::Priors::GenotypeMixturesDirichletAlphaMap& prior_alphas,
                             const LogProbabilityVector& genotype_log_priors,
                             const HaplotypeLikelihoodArray& haplotype_log_likelihoods,
                             const typename SubcloneModelBase<G, GPM>::AlgorithmParameters& params,
                             std::vector<LogProbabilityVector>&& seeds,
                             typename SubcloneModelBase<G, GPM>::OptionalThreadPool workers)
{
    VariationalBayesParameters vb_params {params.epsilon, params.max_iterations};
    if (params.target_max_memory) {
        const auto estimated_memory_default = estimate_memory_requirement<K>(samples, haplotype_log_likelihoods, genotypes.size(), vb_params);
        if (estimated_memory_default > *params.target_max_memory) {
            vb_params.save_memory = true;
        }
    }
    const auto vb_prior_alphas = flatten<K, G, GPM>(prior_alphas, samples);
    const auto log_likelihoods = flatten<K>(genotypes, samples, haplotype_log_likelihoods);
    auto vb_results = octopus::model::run_variational_bayes(vb_prior_alphas, genotype_log_priors, log_likelihoods, vb_params, std::move(seeds), workers);
    return expand<K, G, GPM>(vb_results, samples, std::move(genotype_log_priors));
}

template <typename G, typename GPM,
          std::size_t... Is>
typename SubcloneModelBase<G, GPM>::InferredLatents
run_variational_bayes_helper(const std::vector<SampleName>& samples,
                             const std::vector<G>& genotypes,
                             const typename SubcloneModelBase<G, GPM>::Priors::GenotypeMixturesDirichletAlphaMap& prior_alphas,
                             LogProbabilityVector genotype_log_priors,
                             const HaplotypeLikelihoodArray& haplotype_log_likelihoods,
                             const typename SubcloneModelBase<G, GPM>::AlgorithmParameters& params,
                             std::vector<LogProbabilityVector>&& seeds,
                             typename SubcloneModelBase<G, GPM>::OptionalThreadPool workers,
                             std::index_sequence<Is...>)
{
    constexpr auto max_ploidy = std::index_sequence<Is...>::size();
    const auto ploidy = genotypes.front().ploidy();
    if (ploidy > max_ploidy) {
        throw UnimplementedFeatureError {"ploidies above " + std::to_string(max_ploidy), "SubcloneModel"};
    }
    typename SubcloneModelBase<G, GPM>::InferredLatents result;
    int unused[] = {(ploidy == Is ?
                    (result = run_variational_bayes_helper<Is, G, GPM>(samples, genotypes, prior_alphas, std::move(genotype_log_priors),
                                                                       haplotype_log_likelihoods, params, std::move(seeds), workers)
                                                                       , 0) : 0)...};
    (void) unused;
    return result;
}

template<std::size_t N, std::size_t... Seq>
constexpr std::index_sequence<N + Seq ...>
add(std::index_sequence<Seq...>) { return {}; }

template<std::size_t Min, std::size_t Max>
using make_index_range = decltype(add<Min>(std::make_index_sequence<Max - Min>()));

template <typename G, typename GPM>
typename SubcloneModelBase<G, GPM>::InferredLatents
run_variational_bayes(const std::vector<SampleName>& samples,
                      const std::vector<G>& genotypes,
                      const typename SubcloneModelBase<G, GPM>::Priors& priors,
                      const HaplotypeLikelihoodArray& haplotype_log_likelihoods,
                      const typename SubcloneModelBase<G, GPM>::AlgorithmParameters& params,
                      std::vector<typename SubcloneModelBase<G, GPM>::Latents::LogProbabilityVector> hints,
                      const MappableBlock<Haplotype>* haplotypes,
                      typename SubcloneModelBase<G, GPM>::OptionalThreadPool workers)
{
    constexpr auto max_ploidy = SubcloneModelBase<G, GPM>::max_ploidy;
    auto genotype_log_priors = evaluate(genotypes, priors.genotype_prior_model);
    auto seeds = generate_seeds(samples, genotypes, genotype_log_priors, haplotype_log_likelihoods, priors, params.max_seeds, std::move(hints), haplotypes);
    return run_variational_bayes_helper<G, GPM>(samples, genotypes, priors.alphas, std::move(genotype_log_priors),
                                                haplotype_log_likelihoods, params, std::move(seeds), workers,
                                                make_index_range<1, max_ploidy + 1> {});
}

} // namespace detail

template <typename G, typename GPM>
typename SubcloneModelBase<G, GPM>::InferredLatents
SubcloneModelBase<G, GPM>::evaluate(const std::vector<G>& genotypes,
                                    const HaplotypeLikelihoodArray& haplotype_likelihoods,
                                    std::vector<typename Latents::LogProbabilityVector> hints,
                                    OptionalThreadPool workers) const
{
    assert(!genotypes.empty());
    return detail::run_variational_bayes<G, GPM>(samples_, genotypes, priors_, haplotype_likelihoods, parameters_, std::move(hints), haplotypes_, workers);
}

} // namespace model
} // namespace octopus

#endif
