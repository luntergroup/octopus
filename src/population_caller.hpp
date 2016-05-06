//
//  population_caller.hpp
//  Octopus
//
//  Created by Daniel Cooke on 15/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__population_caller__
#define __Octopus__population_caller__

#include <vector>
#include <string>
#include <memory>

#include "common.hpp"
#include "variant_caller.hpp"
#include "population_genotype_model.hpp"
#include "variant_call.hpp"

class GenomicRegion;
class ReadPipe;
class Variant;
class HaplotypeLikelihoodCache;

namespace Octopus
{
class PopulationVariantCaller : public VariantCaller
{
public:
    struct CallerParameters
    {
        CallerParameters() = default;
        explicit CallerParameters(double min_variant_posterior, double min_refcall_posterior,
                                  unsigned ploidy);
        ~CallerParameters() = default;
        
        double min_variant_posterior;
        double min_refcall_posterior;
        unsigned ploidy;
    };
    
    PopulationVariantCaller() = delete;
    
    explicit PopulationVariantCaller(const ReferenceGenome& reference,
                                     ReadPipe& read_pipe,
                                     CandidateVariantGenerator&& candidate_generator,
                                     VariantCaller::CallerParameters general_parameters,
                                     CallerParameters specific_parameters);
    
    ~PopulationVariantCaller() = default;
    
    PopulationVariantCaller(const PopulationVariantCaller&)            = delete;
    PopulationVariantCaller& operator=(const PopulationVariantCaller&) = delete;
    PopulationVariantCaller(PopulationVariantCaller&&)                 = delete;
    PopulationVariantCaller& operator=(PopulationVariantCaller&&)      = delete;
    
private:
    class Latents : public CallerLatents
    {
    public:
        using ModelInferences = GenotypeModel::Population::InferredLatents;
        
        using CallerLatents::HaplotypeProbabilityMap;
        using CallerLatents::GenotypeProbabilityMap;
        
        friend PopulationVariantCaller;
        
        explicit Latents(const std::vector<SampleIdType>& samples,
                         const std::vector<Haplotype>&,
                         std::vector<Genotype<Haplotype>>&& genotypes,
                         ModelInferences&&);
        explicit Latents(const std::vector<SampleIdType>& samples,
                         const std::vector<Haplotype>&,
                         std::vector<Genotype<Haplotype>>&& genotypes,
                         ModelInferences&&, ModelInferences&&);
        
        std::shared_ptr<HaplotypeProbabilityMap> get_haplotype_posteriors() const noexcept override;
        std::shared_ptr<GenotypeProbabilityMap> get_genotype_posteriors() const noexcept override;
        
    private:
        std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors_;
        std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors_;
        
        boost::optional<ModelInferences> dummy_latents_;
        
        double model_log_evidence_;
        
        HaplotypeProbabilityMap
        calculate_haplotype_posteriors(const std::vector<Haplotype>& haplotypes);
    };
    
    unsigned ploidy_;
    double min_variant_posterior_;
    double min_refcall_posterior_;
    
    std::unique_ptr<CallerLatents>
    infer_latents(const std::vector<Haplotype>& haplotypes,
                  const HaplotypeLikelihoodCache& haplotype_likelihoods) const override;
    
    std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, CallerLatents& latents) const override;
    
    std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, const Latents& latents) const;
    
    std::vector<std::unique_ptr<ReferenceCall>>
    call_reference(const std::vector<Allele>& alleles, CallerLatents& latents,
                   const ReadMap& reads) const override;
};
} // namespace Octopus

#endif /* defined(__Octopus__population_caller__) */
