// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef population_caller_hpp
#define population_caller_hpp

#include <vector>
#include <string>
#include <memory>
#include <typeindex>
#include <map>

#include "config/common.hpp"
#include "basics/ploidy_map.hpp"
#include "basics/phred.hpp"
#include "core/models/mutation/coalescent_model.hpp"
#include "core/models/genotype/genotype_prior_model.hpp"
#include "core/models/genotype/population_prior_model.hpp"
#include "core/models/genotype/independent_population_model.hpp"
#include "core/models/genotype/population_model.hpp"
#include "caller.hpp"

namespace octopus {

class GenomicRegion;
class ReadPipe;
class Variant;
class HaplotypeLikelihoodArray;
class VariantCall;

class PopulationCaller : public Caller
{
public:
    using Caller::CallTypeSet;
    
    struct Parameters
    {
        Phred<double> min_variant_posterior;
        std::vector<unsigned> ploidies;
        boost::optional<CoalescentModel::Parameters> prior_model_params;
        boost::optional<std::size_t> max_genotype_combinations;
        bool use_independent_genotype_priors = false;
        bool deduplicate_haplotypes_with_germline_model = true;
    };
    
    PopulationCaller() = delete;
    
    PopulationCaller(Caller::Components&& components,
                     Caller::Parameters general_parameters,
                     Parameters specific_parameters);
    
    PopulationCaller(const PopulationCaller&)            = delete;
    PopulationCaller& operator=(const PopulationCaller&) = delete;
    PopulationCaller(PopulationCaller&&)                 = delete;
    PopulationCaller& operator=(PopulationCaller&&)      = delete;
    
    ~PopulationCaller() = default;
    
private:
    class Latents;

    using IndexedHaplotypeBlock = MappableBlock<IndexedHaplotype<>>;
    using GenotypeBlock = MappableBlock<Genotype<IndexedHaplotype<>>>;
    
    Parameters parameters_;
    std::vector<unsigned> unique_ploidies_;
    
    std::string do_name() const override;
    CallTypeSet do_call_types() const override;
    unsigned do_min_callable_ploidy() const override;
    unsigned do_max_callable_ploidy() const override;
    
    std::size_t do_remove_duplicates(HaplotypeBlock& haplotypes) const override;
    
    std::unique_ptr<Caller::Latents>
    infer_latents(const HaplotypeBlock& haplotypes,
                  const HaplotypeLikelihoodArray& haplotype_likelihoods,
                  OptionalThreadPool workers) const override;
    
    boost::optional<ModelPosterior>
    calculate_model_posterior(const HaplotypeBlock& haplotypes,
                              const HaplotypeLikelihoodArray& haplotype_likelihoods,
                              const Caller::Latents& latents) const override;
    
    boost::optional<ModelPosterior>
    calculate_model_posterior(const HaplotypeBlock& haplotypes,
                              const HaplotypeLikelihoodArray& haplotype_likelihoods,
                              const Latents& latents) const;
    
    std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, const Caller::Latents& latents, OptionalThreadPool workers) const override;
    
    std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, const Latents& latents) const;
    
    std::vector<std::unique_ptr<ReferenceCall>>
    call_reference(const std::vector<Allele>& alleles, const Caller::Latents& latents,
                   const ReadPileupMap& pileups) const override;
    
    bool use_independence_model() const noexcept;
    std::unique_ptr<Caller::Latents>
    infer_latents_with_joint_model(const HaplotypeBlock& haplotypes,
                                   const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    std::unique_ptr<Caller::Latents>
    infer_latents_with_independence_model(const HaplotypeBlock& haplotypes,
                                          const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    std::unique_ptr<PopulationPriorModel> make_joint_prior_model(const HaplotypeBlock& haplotypes) const;
    std::unique_ptr<GenotypePriorModel> make_independent_prior_model(const HaplotypeBlock& haplotypes) const;

    std::pair<GenotypeBlock, GenotypeBlock>
    propose_model_check_genotypes(std::size_t sample_idx,
                                  const HaplotypeBlock& haplotypes,
                                  const IndexedHaplotypeBlock& indexed_haplotypes,
                                  const Latents& latents) const;
};

class PopulationCaller::Latents : public Caller::Latents
{
public:
    using IndependenceModelInferences = model::IndependentPopulationModel::InferredLatents;
    using ModelInferences = model::PopulationModel::InferredLatents;
    
    using Caller::Latents::HaplotypeProbabilityMap;
    using Caller::Latents::GenotypeProbabilityMap;
    
    friend PopulationCaller;
    
    using IndexedHaplotypeBlock = MappableBlock<IndexedHaplotype<>>;
    using GenotypeBlock = MappableBlock<Genotype<IndexedHaplotype<>>>;
    
    Latents(const std::vector<SampleName>& samples,
            const IndexedHaplotypeBlock&,
            GenotypeBlock genotypes,
            IndependenceModelInferences&&);
    Latents(const std::vector<SampleName>& samples,
            const IndexedHaplotypeBlock&,
            GenotypeBlock,
            ModelInferences&&);
    
    std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors() const noexcept override;
    std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors() const noexcept override;

private:
    GenotypeBlock genotypes_;
    ModelInferences model_latents_;
    std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors_;
    std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors_;
    boost::optional<ModelInferences> dummy_latents_;
};

} // namespace octopus

#endif
