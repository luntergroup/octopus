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

class GenomicRegion;
class ReadPipe;
class Variant;
class VcfRecord;

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
        
        explicit Latents(std::vector<Genotype<Haplotype>>&& germline_genotypes,
                         std::vector<CancerGenotype<Haplotype>>&& somatic_genotypes,
                         GermlineModel::InferredLatents&&,
                         CNVModel::InferredLatents&&,
                         SomaticModel::InferredLatents&&);
        
        std::shared_ptr<HaplotypeProbabilityMap> get_haplotype_posteriors() const override;
        std::shared_ptr<GenotypeProbabilityMap> get_genotype_posteriors() const override;
        
    private:
        std::vector<Genotype<Haplotype>> germline_genotypes_;
        std::vector<CancerGenotype<Haplotype>> somatic_genotypes_;
        
        GermlineModel::InferredLatents germline_model_inferences_;
        CNVModel::InferredLatents cnv_model_inferences_;
        SomaticModel::InferredLatents somatic_model_inferences_;
        
        friend CancerVariantCaller;
    };
    
    struct ModelPosteriors
    {
        explicit ModelPosteriors(double germline, double cnv, double somatic);
        double germline, cnv, somatic;
    };
    
    CallerParameters parameters_;
    
    std::unique_ptr<CallerLatents>
    infer_latents(const std::vector<Haplotype>& haplotypes,
                  const HaplotypeLikelihoodCache& haplotype_likelihoods) const override;
    
    CNVModel::Priors calculate_cnv_model_priors(const CoalescentModel& prior_model) const;
    SomaticModel::Priors calculate_somatic_model_priors(const SomaticMutationModel& prior_model) const;
    
    std::vector<VcfRecord::Builder>
    call_variants(const std::vector<Variant>& candidates,
                  const std::vector<Allele>& callable_alleles,
                  CallerLatents* latents,
                  const Phaser::PhaseSet& phase_set,
                  const ReadMap& reads) const override;
    
    ModelPosteriors calculate_model_posteriors(const Latents& inferences) const;
    
    std::vector<VcfRecord::Builder>
    call_germline_variants(const std::vector<Variant>& candidates,
                           const std::vector<Allele>& callable_alleles,
                           const GermlineModel::Latents& posteriors,
                           const Phaser::PhaseSet& phase_set,
                           const ReadMap& reads) const;
    std::vector<VcfRecord::Builder>
    call_cnv_variants(const std::vector<Variant>& candidates,
                      const std::vector<Allele>& callable_alleles,
                      const CNVModel::Latents& posteriors,
                      const Phaser::PhaseSet& phase_set,
                      const ReadMap& reads) const;
    std::vector<VcfRecord::Builder>
    call_somatic_variants(const std::vector<Variant>& candidates,
                          const std::vector<Allele>& callable_alleles,
                          const SomaticModel::Latents& posteriors,
                          const Phaser::PhaseSet& phase_set,
                          const ReadMap& reads) const;
};

} // namespace Octopus

#endif /* defined(__Octopus__cancer_caller__) */
