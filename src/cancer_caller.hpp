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

#include <boost/optional.hpp>

#include "common.hpp"
#include "variant_caller.hpp"
#include "coalescent_model.hpp"
#include "somatic_mutation_model.hpp"
#include "individual_genotype_model.hpp"
#include "cnv_genotype_model.hpp"
#include "somatic_genotype_model.hpp"
#include "variant_call.hpp"

class GenomicRegion;
class ReadPipe;
class Variant;

namespace Octopus
{
class CancerVariantCaller : public VariantCaller
{
public:
    struct CallerParameters
    {
        CallerParameters() = default;
        explicit CallerParameters(double min_variant_posterior, double min_somatic_posterior,
                                  double min_refcall_posterior, unsigned ploidy,
                                  boost::optional<SampleIdType> normal_sample,
                                  double somatic_mutation_rate, bool call_somatics_only);
        ~CallerParameters() = default;
        
        double min_variant_posterior;
        double min_somatic_posterior;
        double min_refcall_posterior;
        unsigned ploidy;
        boost::optional<SampleIdType> normal_sample;
        double somatic_mutation_rate;
        bool call_somatics_only;
        std::size_t max_genotypes;
    };
    
    CancerVariantCaller() = delete;
    
    explicit CancerVariantCaller(const ReferenceGenome& reference,
                                 ReadPipe& read_pipe,
                                 CandidateVariantGenerator&& candidate_generator,
                                 VariantCaller::CallerParameters general_parameters,
                                 CallerParameters specific_parameters);
    
    ~CancerVariantCaller() = default;
    
    CancerVariantCaller(const CancerVariantCaller&)            = delete;
    CancerVariantCaller& operator=(const CancerVariantCaller&) = delete;
    CancerVariantCaller(CancerVariantCaller&&)                 = delete;
    CancerVariantCaller& operator=(CancerVariantCaller&&)      = delete;
    
private:
    using GermlineModel = GenotypeModel::Individual;
    using CNVModel      = GenotypeModel::CNV;
    using SomaticModel  = GenotypeModel::Somatic;
    
    class Latents : public CallerLatents
    {
    public:
        using CallerLatents::HaplotypeProbabilityMap;
        using CallerLatents::GenotypeProbabilityMap;
        
        struct NormalInferences
        {
            NormalInferences() = default;
            GermlineModel::InferredLatents germline, dummy;
        };
        
        explicit Latents(std::vector<Genotype<Haplotype>>&& germline_genotypes,
                         std::vector<CancerGenotype<Haplotype>>&& somatic_genotypes,
                         GermlineModel::InferredLatents&&,
                         CNVModel::InferredLatents&&,
                         SomaticModel::InferredLatents&&);
        
        explicit Latents(std::vector<Genotype<Haplotype>>&& germline_genotypes,
                         std::vector<CancerGenotype<Haplotype>>&& somatic_genotypes,
                         GermlineModel::InferredLatents&&,
                         CNVModel::InferredLatents&&,
                         SomaticModel::InferredLatents&&,
                         NormalInferences&&);
        
        std::shared_ptr<HaplotypeProbabilityMap> get_haplotype_posteriors() const override;
        std::shared_ptr<GenotypeProbabilityMap> get_genotype_posteriors() const override;
        
    private:
        std::vector<Genotype<Haplotype>> germline_genotypes_;
        std::vector<CancerGenotype<Haplotype>> somatic_genotypes_;
        
        GermlineModel::InferredLatents germline_model_inferences_;
        CNVModel::InferredLatents cnv_model_inferences_;
        SomaticModel::InferredLatents somatic_model_inferences_;
        
        boost::optional<NormalInferences> normal_validation_inferences_;
        
        friend CancerVariantCaller;
    };
    
    struct ModelPosteriors
    {
        explicit ModelPosteriors(double germline, double cnv, double somatic);
        double germline, cnv, somatic;
    };
    
    CallerParameters parameters_;
    
    // overrides
    
    std::unique_ptr<CallerLatents>
    infer_latents(const std::vector<Haplotype>& haplotypes,
                  const HaplotypeLikelihoodCache& haplotype_likelihoods) const override;
    
    std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, CallerLatents& latents) const override;
    
    std::vector<std::unique_ptr<ReferenceCall>>
    call_reference(const std::vector<Allele>& alleles, CallerLatents& latents,
                   const ReadMap& reads) const override;
    
    // other private methods
    
    using GermlineGenotypeReference      = Genotype<Haplotype>;
    using GermlineGenotypeProbabilityMap = std::unordered_map<GermlineGenotypeReference, double>;
    using ProbabilityVector = std::vector<double>;
    
    bool has_normal_sample() const noexcept;
    const SampleIdType& normal_sample() const;
    
    void filter(std::vector<CancerGenotype<Haplotype>>& cancer_genotypes,
                const std::vector<Genotype<Haplotype>>& germline_genotypes,
                const GermlineModel::InferredLatents& germline_inferences,
                const CNVModel::InferredLatents& cnv_inferences) const;
    
    std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, Latents& latents) const;
    
    CNVModel::Priors calculate_cnv_model_priors(const CoalescentModel& prior_model) const;
    SomaticModel::Priors calculate_somatic_model_priors(const SomaticMutationModel& prior_model) const;

    ModelPosteriors calculate_model_posteriors(const Latents& inferences) const;
    
    GermlineGenotypeProbabilityMap
    calculate_germline_genotype_posteriors(const Latents& inferences,
                                           const ModelPosteriors& model_posteriors) const;
    
    ProbabilityVector
    calculate_probability_samples_not_somatic(const Latents& inferences) const;
    
    double calculate_somatic_probability(const ProbabilityVector& sample_somatic_posteriors,
                                         const ModelPosteriors& model_posteriors) const;
};

} // namespace Octopus

#endif /* defined(__Octopus__cancer_caller__) */
