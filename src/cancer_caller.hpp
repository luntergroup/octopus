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

#include "common.hpp"
#include "variant_caller.hpp"
#include "cancer_genotype_model.hpp"

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
                                  SampleIdType normal_sample, bool call_somatics_only);
        ~CallerParameters() = default;
        
        double min_variant_posterior;
        double min_somatic_posterior;
        double min_refcall_posterior;
        unsigned ploidy;
        SampleIdType normal_sample;
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
    class Latents : public CallerLatents
    {
    public:
        using ModelLatents = GenotypeModel::Cancer::InferredLatents;
        
        using CallerLatents::HaplotypePosteriorMap;
        using CallerLatents::GenotypePosteriorMap;
        
        explicit Latents(ModelLatents&&);
        
        std::shared_ptr<HaplotypePosteriorMap> get_haplotype_posteriors() const override;
        std::shared_ptr<GenotypePosteriorMap> get_genotype_posteriors() const override;
        
    private:
        ModelLatents model_latents_;
    };
    
    CallerParameters parameters_;
    
    std::unique_ptr<CallerLatents>
    infer_latents(const std::vector<Haplotype>& haplotypes,
                  const HaplotypeLikelihoodCache& haplotype_likelihoods) const override;
    
    std::vector<VcfRecord::Builder>
    call_variants(const std::vector<Variant>& candidates,
                  const std::vector<Allele>& callable_alleles,
                  CallerLatents* latents,
                  const Phaser::PhaseSet& phase_set,
                  const ReadMap& reads) const override;
};

} // namespace Octopus

#endif /* defined(__Octopus__cancer_caller__) */
