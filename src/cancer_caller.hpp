//
//  cancer_caller.hpp
//  Octopus
//
//  Created by Daniel Cooke on 16/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__cancer_caller__
#define __Octopus__cancer_caller__

#include <vector>
#include <unordered_map>
#include <memory>
#include <functional>
#include <typeindex>

#include <boost/optional.hpp>

#include "common.hpp"
#include "variant_caller.hpp"
#include "coalescent_model.hpp"
#include "somatic_mutation_model.hpp"
#include "individual_genotype_model.hpp"
#include "cnv_genotype_model.hpp"
#include "somatic_genotype_model.hpp"
#include "variant_call.hpp"
#include "phred.hpp"

class GenomicRegion;
class ReadPipe;
class Variant;

namespace octopus
{
class CancerVariantCaller : public VariantCaller
{
public:
    using VariantCaller::CallTypeSet;
    using VariantCaller::CallerComponents;
    
    struct CallerParameters
    {
        Phred<double> min_variant_posterior;
        Phred<double> min_somatic_posterior;
        Phred<double> min_refcall_posterior;
        unsigned ploidy;
        boost::optional<SampleName> normal_sample;
        double somatic_mutation_rate;
        double min_somatic_frequency;
        double credible_mass;
        unsigned max_genotypes;
    };
    
    CancerVariantCaller() = delete;
    
    CancerVariantCaller(CallerComponents&& components,
                        VariantCaller::CallerParameters general_parameters,
                        CallerParameters specific_parameters);
    
    CancerVariantCaller(const CancerVariantCaller&)            = delete;
    CancerVariantCaller& operator=(const CancerVariantCaller&) = delete;
    CancerVariantCaller(CancerVariantCaller&&)                 = delete;
    CancerVariantCaller& operator=(CancerVariantCaller&&)      = delete;
    
    ~CancerVariantCaller() = default;
    
private:
    using GermlineModel = model::Individual;
    using CNVModel      = model::CNV;
    using SomaticModel  = model::Somatic;
    
    class Latents;
    
    friend Latents;
    
    struct ModelProbabilities
    {
        double germline, cnv, somatic;
    };
    
    using ModelPriors     = ModelProbabilities;
    using ModelPosteriors = ModelProbabilities;
    
    CallerParameters parameters_;
    
    // overrides
    
    CallTypeSet do_get_call_types() const override;
    
    std::unique_ptr<CallerLatents>
    infer_latents(const std::vector<Haplotype>& haplotypes,
                  const HaplotypeLikelihoodCache& haplotype_likelihoods) const override;
    
    boost::optional<double>
    calculate_dummy_model_posterior(const std::vector<Haplotype>& haplotypes,
                                    const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                    const CallerLatents& latents) const override;
    
    boost::optional<double>
    calculate_dummy_model_posterior(const std::vector<Haplotype>& haplotypes,
                                    const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                    const Latents& latents) const;
    
    std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, const CallerLatents& latents) const override;
    
    std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, const Latents& latents) const;
    
    std::vector<std::unique_ptr<ReferenceCall>>
    call_reference(const std::vector<Allele>& alleles, const CallerLatents& latents,
                   const ReadMap& reads) const override;
    
    // helpers
    
    bool has_normal_sample() const noexcept;
    const SampleName& normal_sample() const;
    
    void filter(std::vector<CancerGenotype<Haplotype>>& cancer_genotypes,
                const std::vector<Genotype<Haplotype>>& germline_genotypes,
                const GermlineModel::InferredLatents& germline_inferences,
                const CNVModel::InferredLatents& cnv_inferences) const;
    
    using GermlineGenotypeReference      = Genotype<Haplotype>;
    using GermlineGenotypeProbabilityMap = std::unordered_map<GermlineGenotypeReference, double>;
    using ProbabilityVector = std::vector<double>;
    
    CNVModel::Priors get_cnv_model_priors(const CoalescentModel& prior_model) const;
    SomaticModel::Priors get_somatic_model_priors(const SomaticMutationModel& prior_model) const;
    
    ModelPriors get_model_priors() const;
    
    ModelPosteriors calculate_model_posteriors(const Latents& inferences) const;
    
    GermlineGenotypeProbabilityMap calculate_germline_genotype_posteriors(const Latents& inferences,
                                                                          const ModelPosteriors& model_posteriors) const;
    
    ProbabilityVector calculate_probability_samples_not_somatic(const Latents& inferences) const;
    
    Phred<double> calculate_somatic_probability(const ProbabilityVector& sample_somatic_posteriors,
                                                const ModelPosteriors& model_posteriors) const;
};

class CancerVariantCaller::Latents : public CallerLatents
{
public:
    using CallerLatents::HaplotypeProbabilityMap;
    using CallerLatents::GenotypeProbabilityMap;
    
    Latents() = delete;
    
    Latents(const std::vector<Haplotype>& haplotypes,
            CancerVariantCaller::ModelPriors model_priors,
            std::vector<Genotype<Haplotype>>&& germline_genotypes,
            std::vector<CancerGenotype<Haplotype>>&& somatic_genotypes,
            GermlineModel::InferredLatents&&, CNVModel::InferredLatents&&,
            SomaticModel::InferredLatents&&,
            const std::vector<SampleName>& samples,
            boost::optional<std::reference_wrapper<const SampleName>> normal_sample);
    
    std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors() const override;
    std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors() const override;
    
private:
    std::vector<Genotype<Haplotype>> germline_genotypes_;
    std::vector<CancerGenotype<Haplotype>> somatic_genotypes_;
    
    CancerVariantCaller::ModelPriors model_priors_;
    
    GermlineModel::InferredLatents germline_model_inferences_;
    CNVModel::InferredLatents cnv_model_inferences_;
    SomaticModel::InferredLatents somatic_model_inferences_;
    
    std::reference_wrapper<const std::vector<Haplotype>> haplotypes_;
    
    std::reference_wrapper<const std::vector<SampleName>> samples_;
    boost::optional<std::reference_wrapper<const SampleName>> normal_sample_;
    
    friend CancerVariantCaller;
};

} // namespace octopus

#endif /* defined(__Octopus__cancer_caller__) */
