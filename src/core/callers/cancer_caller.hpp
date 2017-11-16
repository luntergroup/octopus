// Copyright (c) 2017 Daniel Cooke
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
#include "core/types/haplotype.hpp"
#include "core/types/cancer_genotype.hpp"
#include "core/models/genotype/cancer_genotype_prior_model.hpp"
#include "core/models/genotype/genotype_prior_model.hpp"
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
        boost::optional<CoalescentModel::Parameters> germline_prior_model_params;
        SomaticMutationModel::Parameters somatic_mutation_model_params;
        double min_expected_somatic_frequency, credible_mass, min_credible_somatic_frequency;
        unsigned max_genotypes = 20000;
        double cnv_normal_alpha = 10.0, cnv_tumour_alpha = 0.75;
        double somatic_normal_germline_alpha = 10.0, somatic_normal_somatic_alpha = 0.08;
        double somatic_tumour_germline_alpha = 1.0, somatic_tumour_somatic_alpha = 0.8;
        double germline_weight = 70, cnv_weight = 3, somatic_weight = 2;
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
    
    bool has_normal_sample() const noexcept;
    const SampleName& normal_sample() const;
    
    using GenotypeVector                 = std::vector<Genotype<Haplotype>>;
    using CancerGenotypeVector           = std::vector<CancerGenotype<Haplotype>>;
    using GermlineGenotypeReference      = Genotype<Haplotype>;
    using GermlineGenotypeProbabilityMap = std::unordered_map<GermlineGenotypeReference, double>;
    using ProbabilityVector              = std::vector<double>;
    
    void evaluate_germline_model(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    void evaluate_cnv_model(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    void evaluate_tumour_model(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    void evaluate_noise_model(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    
    std::unique_ptr<GenotypePriorModel> make_germline_prior_model(const std::vector<Haplotype>& haplotypes) const;
    CNVModel::Priors get_cnv_model_priors(const GenotypePriorModel& prior_model) const;
    TumourModel::Priors get_somatic_model_priors(const CancerGenotypePriorModel& prior_model) const;
    TumourModel::Priors get_noise_model_priors(const CancerGenotypePriorModel& prior_model) const;
    CNVModel::Priors get_normal_noise_model_priors(const GenotypePriorModel& prior_model) const;
    
    CancerGenotypeVector generate_cancer_genotypes(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    
    ModelPriors get_model_priors() const;
    ModelPosteriors calculate_model_posteriors(const Latents& latents) const;
    
    GermlineGenotypeProbabilityMap calculate_germline_genotype_posteriors(const Latents& latents,
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
            const std::vector<SampleName>& samples,
            CancerCaller::ModelPriors model_priors);
    
    std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors() const override;
    std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors() const override;
    
private:
    std::reference_wrapper<const std::vector<Haplotype>> haplotypes_;
    std::vector<Genotype<Haplotype>> germline_genotypes_;
    std::vector<CancerGenotype<Haplotype>> cancer_genotypes_;
    
    std::reference_wrapper<const std::vector<SampleName>> samples_;
    boost::optional<std::reference_wrapper<const SampleName>> normal_sample_ = boost::none;
    
    CancerCaller::ModelPriors model_priors_;
    std::unique_ptr<GenotypePriorModel> germline_prior_model_ = nullptr;
    boost::optional<CancerGenotypePriorModel> cancer_genotype_prior_model_ = boost::none;
    std::unique_ptr<GermlineModel> germline_model_ = nullptr;
    GermlineModel::InferredLatents germline_model_inferences_;
    CNVModel::InferredLatents cnv_model_inferences_;
    TumourModel::InferredLatents tumour_model_inferences_;
    boost::optional<TumourModel::InferredLatents> noise_model_inferences_ = boost::none;
    boost::optional<GermlineModel::InferredLatents> normal_germline_inferences_ = boost::none;
    
    mutable std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors_ = nullptr;
    mutable std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors_ = nullptr;
    
    friend CancerCaller;
    
    void compute_genotype_posteriors() const;
    void compute_haplotype_posteriors() const;
};

} // namespace octopus

#endif
