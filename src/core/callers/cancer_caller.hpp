// Copyright (c) 2015-2018 Daniel Cooke
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
#include "core/models/genotype/subclone_model.hpp"
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
        enum class NormalContaminationRisk { high, low };
        Phred<double> min_variant_posterior, min_somatic_posterior, min_refcall_posterior;
        unsigned ploidy;
        boost::optional<SampleName> normal_sample;
        boost::optional<CoalescentModel::Parameters> germline_prior_model_params;
        SomaticMutationModel::Parameters somatic_mutation_model_params;
        double min_expected_somatic_frequency, credible_mass, min_credible_somatic_frequency;
        std::size_t max_genotypes = 20000;
        unsigned max_somatic_haplotypes = 1;
        NormalContaminationRisk normal_contamination_risk = NormalContaminationRisk::low;
        double cnv_normal_alpha = 50.0, cnv_tumour_alpha = 0.5;
        double somatic_normal_germline_alpha = 50.0, somatic_normal_somatic_alpha = 0.05;
        double somatic_tumour_germline_alpha = 1.5, somatic_tumour_somatic_alpha = 1.0;
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
    using CNVModel      = model::SubcloneModel;
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
    
    void generate_germline_genotypes(Latents& latents, const std::vector<Haplotype>& haplotypes) const;
    void generate_cancer_genotypes(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    void generate_cancer_genotypes_with_clean_normal(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    void generate_cancer_genotypes_with_contaminated_normal(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    void generate_cancer_genotypes_with_no_normal(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    void generate_cancer_genotypes(Latents& latents, const std::vector<Genotype<Haplotype>>& germline_genotypes) const;
    bool has_high_normal_contamination_risk(const Latents& latents) const;
    
    void evaluate_germline_model(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    void evaluate_cnv_model(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    void evaluate_tumour_model(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    void evaluate_noise_model(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    
    void set_model_priors(Latents& latents) const;
    void set_model_posteriors(Latents& latents) const;
    
    void fit_tumour_model(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    
    std::unique_ptr<GenotypePriorModel> make_germline_prior_model(const std::vector<Haplotype>& haplotypes) const;
    CNVModel::Priors get_cnv_model_priors(const GenotypePriorModel& prior_model) const;
    TumourModel::Priors get_somatic_model_priors(const CancerGenotypePriorModel& prior_model, unsigned somatic_ploidy) const;
    TumourModel::Priors get_noise_model_priors(const CancerGenotypePriorModel& prior_model, unsigned somatic_ploidy) const;
    CNVModel::Priors get_normal_noise_model_priors(const GenotypePriorModel& prior_model) const;
    
    GermlineGenotypeProbabilityMap calculate_germline_genotype_posteriors(const Latents& latents) const;
    Phred<double> calculate_somatic_probability(const Latents& latents) const;
    
    // logging
    void log(const ModelPosteriors& model_posteriors) const;
    void log(const GenotypeVector& germline_genotypes,
             const GermlineGenotypeProbabilityMap& germline_genotype_posteriors,
             const CNVModel::InferredLatents& cnv_inferences,
             const CancerGenotypeVector& cancer_genotypes,
             const TumourModel::InferredLatents& tumour_inferences) const;
};

class CancerCaller::Latents : public Caller::Latents
{
public:
    using Caller::Latents::HaplotypeProbabilityMap;
    using Caller::Latents::GenotypeProbabilityMap;
    
    Latents() = delete;
    
    Latents(const std::vector<Haplotype>& haplotypes,
            const std::vector<SampleName>& samples);
    
    std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors() const override;
    std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors() const override;
    
private:
    std::reference_wrapper<const std::vector<Haplotype>> haplotypes_;
    std::vector<Genotype<Haplotype>> germline_genotypes_;
    unsigned somatic_ploidy_ = 1;
    std::vector<CancerGenotype<Haplotype>> cancer_genotypes_;
    boost::optional<std::vector<std::vector<unsigned>>> germline_genotype_indices_ = boost::none;
    boost::optional<std::vector<CancerGenotypeIndex>> cancer_genotype_indices_ = boost::none;
    
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
    CancerCaller::ModelPosteriors model_posteriors_;
    
    mutable std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors_ = nullptr;
    mutable std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors_ = nullptr;
    
    friend CancerCaller;
    
    void compute_genotype_posteriors() const;
    void compute_haplotype_posteriors() const;
};

} // namespace octopus

#endif
