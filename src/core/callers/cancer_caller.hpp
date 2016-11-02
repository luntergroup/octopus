// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef cancer_caller_hpp
#define cancer_caller_hpp

#include <vector>
#include <unordered_map>
#include <memory>
#include <functional>
#include <typeindex>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "core/models/mutation/coalescent_model.hpp"
#include "core/models/mutation/somatic_mutation_model.hpp"
#include "core/models/genotype/individual_model.hpp"
#include "core/models/genotype/cnv_model.hpp"
#include "core/models/genotype/tumour_model.hpp"
#include "basics/phred.hpp"
#include "caller.hpp"

namespace octopus {

class GenomicRegion;
class ReadPipe;
class Variant;
class VariantCall;

class CancerCaller : public Caller
{
public:
    using Caller::CallTypeSet;
    
    struct Parameters
    {
        Phred<double> min_variant_posterior, min_somatic_posterior, min_refcall_posterior;
        unsigned ploidy;
        boost::optional<SampleName> normal_sample;
        CoalescentModel::Parameters germline_prior_model_params;
        SomaticMutationModel::Parameters somatic_mutation_model_params;
        double min_somatic_frequency, credible_mass;
        unsigned max_genotypes = 30000;
        double cnv_normal_alpha = 50.0, cnv_tumour_alpha = 0.75;
        double somatic_normal_germline_alpha = 50.0, somatic_normal_somatic_alpha = 0.01;
        double somatic_tumour_germline_alpha = 5.0, somatic_tumour_somatic_alpha = 0.8;
        double germline_weight = 90, cnv_weight = 5, somatic_weight = 1;
    };
    
    CancerCaller() = delete;
    
    CancerCaller(Caller::Components&& components,
                 Caller::Parameters general_parameters,
                 Parameters specific_parameters);
    
    CancerCaller(const CancerCaller&)            = delete;
    CancerCaller& operator=(const CancerCaller&) = delete;
    CancerCaller(CancerCaller&&)                 = delete;
    CancerCaller& operator=(CancerCaller&&)      = delete;
    
    ~CancerCaller() = default;
    
private:
    using GermlineModel = model::IndividualModel;
    using CNVModel      = model::CNVModel;
    using TumourModel   = model::TumourModel;
    
    class Latents;
    
    friend Latents;
    
    struct ModelProbabilities
    {
        double germline, cnv, somatic;
    };
    
    using ModelPriors     = ModelProbabilities;
    using ModelPosteriors = ModelProbabilities;
    
    Parameters parameters_;
    
    // overrides
    
    std::string do_name() const override;
    CallTypeSet do_call_types() const override;
    
    std::unique_ptr<Caller::Latents>
    infer_latents(const std::vector<Haplotype>& haplotypes,
                  const HaplotypeLikelihoodCache& haplotype_likelihoods) const override;
    
    boost::optional<double>
    calculate_model_posterior(const std::vector<Haplotype>& haplotypes,
                              const HaplotypeLikelihoodCache& haplotype_likelihoods,
                              const Caller::Latents& latents) const override;
    
    boost::optional<double>
    calculate_model_posterior(const std::vector<Haplotype>& haplotypes,
                              const HaplotypeLikelihoodCache& haplotype_likelihoods,
                              const Latents& latents) const;
    
    std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, const Caller::Latents& latents) const override;
    
    std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, const Latents& latents) const;
    
    std::vector<std::unique_ptr<ReferenceCall>>
    call_reference(const std::vector<Allele>& alleles, const Caller::Latents& latents,
                   const ReadMap& reads) const override;
    
    // helpers
    
    bool has_normal_sample() const noexcept;
    const SampleName& normal_sample() const;
    
    void reduce(std::vector<CancerGenotype<Haplotype>>& cancer_genotypes,
                const std::vector<Genotype<Haplotype>>& germline_genotypes,
                const GermlineModel& germline_model,
                const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    
    using GermlineGenotypeReference      = Genotype<Haplotype>;
    using GermlineGenotypeProbabilityMap = std::unordered_map<GermlineGenotypeReference, double>;
    using ProbabilityVector = std::vector<double>;
    
    CNVModel::Priors get_cnv_model_priors(const CoalescentModel& prior_model) const;
    TumourModel::Priors get_somatic_model_priors(const SomaticMutationModel& prior_model) const;
    ModelPriors get_model_priors() const;
    ModelPosteriors calculate_model_posteriors(const Latents& inferences) const;
    GermlineGenotypeProbabilityMap calculate_germline_genotype_posteriors(const Latents& inferences,
                                                                          const ModelPosteriors& model_posteriors) const;
    ProbabilityVector calculate_probability_samples_not_somatic(const Latents& inferences) const;
    Phred<double> calculate_somatic_probability(const ProbabilityVector& sample_somatic_posteriors,
                                                const ModelPosteriors& model_posteriors) const;
};

class CancerCaller::Latents : public Caller::Latents
{
public:
    using Caller::Latents::HaplotypeProbabilityMap;
    using Caller::Latents::GenotypeProbabilityMap;
    
    Latents() = delete;
    
    Latents(const std::vector<Haplotype>& haplotypes,
            CancerCaller::ModelPriors model_priors,
            std::vector<Genotype<Haplotype>>&& germline_genotypes,
            std::vector<CancerGenotype<Haplotype>>&& somatic_genotypes,
            GermlineModel::InferredLatents&&, CNVModel::InferredLatents&&,
            TumourModel::InferredLatents&&,
            const std::vector<SampleName>& samples,
            boost::optional<std::reference_wrapper<const SampleName>> normal_sample);
    
    std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors() const override;
    std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors() const override;
    
private:
    std::vector<Genotype<Haplotype>> germline_genotypes_;
    std::vector<CancerGenotype<Haplotype>> somatic_genotypes_;
    
    CancerCaller::ModelPriors model_priors_;
    
    GermlineModel::InferredLatents germline_model_inferences_;
    CNVModel::InferredLatents cnv_model_inferences_;
    TumourModel::InferredLatents somatic_model_inferences_;
    
    std::reference_wrapper<const std::vector<Haplotype>> haplotypes_;
    
    std::reference_wrapper<const std::vector<SampleName>> samples_;
    boost::optional<std::reference_wrapper<const SampleName>> normal_sample_;
    
    friend CancerCaller;
};

} // namespace octopus

#endif
