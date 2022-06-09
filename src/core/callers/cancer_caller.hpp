// Copyright (c) 2015-2021 Daniel Cooke
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
        
        struct Concentrations
        {
            struct CNVModelConcentrations
            {
                double normal = 50., tumour = .5;
            };
            struct SomaticModelConcentrations
            {
                double normal_germline = 50., normal_somatic = .05, tumour_germline = 1.5, tumour_somatic = 1.;
            };
            
            CNVModelConcentrations cnv;
            SomaticModelConcentrations somatic;
        };
        
        Phred<double> min_variant_posterior, min_somatic_posterior;
        unsigned ploidy;
        boost::optional<SampleName> normal_sample;
        boost::optional<CoalescentModel::Parameters> germline_prior_model_params;
        SomaticMutationModel::Parameters somatic_mutation_model_params;
        double min_expected_somatic_frequency, credible_mass, min_credible_somatic_frequency;
        boost::optional<std::size_t> max_genotypes = 20000;
        unsigned max_somatic_haplotypes = 1;
        NormalContaminationRisk normal_contamination_risk = NormalContaminationRisk::low;
        bool deduplicate_haplotypes_with_germline_model = true;
        boost::optional<unsigned> max_vb_seeds = boost::none; // Use default if none
        Concentrations concentrations = Concentrations {};
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
    using SomaticModel  = model::SomaticSubcloneModel;
    
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
    
    bool has_normal_sample() const noexcept;
    const SampleName& normal_sample() const;
    
    using IndexedHaplotypeBlock          = MappableBlock<IndexedHaplotype<>>;
    using GenotypeVector                 = MappableBlock<Genotype<IndexedHaplotype<>>>;
    using CancerGenotypeVector           = MappableBlock<CancerGenotype<IndexedHaplotype<>>>;
    using GermlineGenotypeReference      = Genotype<IndexedHaplotype<>>;
    using GermlineGenotypeProbabilityMap = std::unordered_map<GermlineGenotypeReference, double>;
    using ProbabilityVector              = std::vector<double>;
    
    void generate_germline_genotypes(Latents& latents, const IndexedHaplotypeBlock& haplotypes) const;
    void generate_cancer_genotypes(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    void generate_cancer_genotypes_with_clean_normal(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    void generate_cancer_genotypes_with_contaminated_normal(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    void generate_cancer_genotypes_with_no_normal(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    void generate_cancer_genotypes(Latents& latents, const MappableBlock<Genotype<IndexedHaplotype<>>>& germline_genotypes) const;
    bool has_high_normal_contamination_risk(const Latents& latents) const;
    
    void evaluate_germline_model(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    void evaluate_cnv_model(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods, OptionalThreadPool workers) const;
    void evaluate_somatic_model(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods, OptionalThreadPool workers) const;
    void evaluate_noise_model(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    
    void set_model_priors(Latents& latents) const;
    void set_model_posteriors(Latents& latents) const;

    void set_cancer_genotype_prior_model(Latents& latents) const;
    void fit_somatic_model(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods, OptionalThreadPool workers) const;
    
    std::unique_ptr<GenotypePriorModel> make_germline_prior_model(const HaplotypeBlock& haplotypes) const;
    CNVModel::Priors get_cnv_model_priors(const GenotypePriorModel& prior_model) const;
    SomaticModel::Priors get_somatic_model_priors(const CancerGenotypePriorModel& prior_model, unsigned somatic_ploidy) const;
    SomaticModel::Priors get_noise_model_priors(const CancerGenotypePriorModel& prior_model, unsigned somatic_ploidy) const;
    CNVModel::Priors get_normal_noise_model_priors(const GenotypePriorModel& prior_model) const;
    
    GermlineGenotypeProbabilityMap calculate_germline_genotype_posteriors(const Latents& latents) const;
    double calculate_somatic_mass(const Latents& latents) const;
    Phred<double> calculate_segregation_probability(const Variant& variant, const Latents& latents, double somatic_mass) const;
    
    // logging
    void log(const ModelPosteriors& model_posteriors) const;
    void log(const GenotypeVector& germline_genotypes,
             const GermlineGenotypeProbabilityMap& germline_genotype_posteriors,
             const GermlineModel::InferredLatents& germline_inferences,
             const CNVModel::InferredLatents& cnv_inferences,
             const CancerGenotypeVector& cancer_genotypes,
             const SomaticModel::InferredLatents& tumour_inferences) const;
};

class CancerCaller::Latents : public Caller::Latents
{
public:
    using Caller::Latents::HaplotypeProbabilityMap;
    using Caller::Latents::GenotypeProbabilityMap;
    
    Latents() = delete;
    
    Latents(const HaplotypeBlock& haplotypes,
            const std::vector<SampleName>& samples,
            const CancerCaller::Parameters& parameters);
    
    std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors() const override;
    std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors() const override;
    
private:
    std::reference_wrapper<const HaplotypeBlock> haplotypes_;
    MappableBlock<IndexedHaplotype<>> indexed_haplotypes_;
    std::reference_wrapper<const std::vector<SampleName>> samples_;
    std::reference_wrapper<const CancerCaller::Parameters> parameters_;
    
    // prior stuff
    CancerCaller::ModelPriors model_priors_;
    std::unique_ptr<GenotypePriorModel> germline_prior_model_ = nullptr;
    boost::optional<CancerGenotypePriorModel> cancer_genotype_prior_model_ = boost::none;
    std::unique_ptr<GermlineModel> germline_model_ = nullptr;
    // germline and CNV model
    MappableBlock<Genotype<IndexedHaplotype<>>> germline_genotypes_;
    GermlineModel::InferredLatents germline_model_inferences_;
    CNVModel::InferredLatents cnv_model_inferences_;
    // somatic model
    std::vector<MappableBlock<CancerGenotype<IndexedHaplotype<>>>> cancer_genotypes_;
    std::vector<SomaticModel::InferredLatents> somatic_model_inferences_;
    std::vector<double> somatic_model_posteriors_;
    std::size_t max_evidence_somatic_model_index_;
    unsigned inferred_somatic_ploidy_ = 1;
    // noise model
    boost::optional<SomaticModel::InferredLatents> noise_model_inferences_ = boost::none;
    boost::optional<GermlineModel::InferredLatents> normal_germline_inferences_ = boost::none;
    // posterior stuff
    CancerCaller::ModelPosteriors model_posteriors_;
    
    mutable std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors_ = nullptr;
    mutable std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors_ = nullptr;
    
    friend CancerCaller;
    
    void compute_genotype_posteriors() const;
    void compute_haplotype_posteriors() const;
};

} // namespace octopus

#endif
