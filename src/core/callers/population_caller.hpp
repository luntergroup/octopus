// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef __Octopus__population_caller__
#define __Octopus__population_caller__

#include <vector>
#include <string>
#include <memory>
#include <typeindex>

#include <config/common.hpp>
#include <core/models/genotype/population_model.hpp>
#include <basics/phred.hpp>

#include "caller.hpp"

namespace octopus {

class GenomicRegion;
class ReadPipe;
class Variant;
class HaplotypeLikelihoodCache;
class VariantCall;

class PopulationCaller : public Caller
{
public:
    using Caller::CallTypeSet;
    
    struct Parameters
    {
        Phred<double> min_variant_posterior;
        Phred<double> min_refcall_posterior;
        unsigned ploidy;
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
    class Latents : public Caller::Latents
    {
    public:
        using ModelInferences = model::PopulationModel::InferredLatents;
        
        using Caller::Latents::HaplotypeProbabilityMap;
        using Caller::Latents::GenotypeProbabilityMap;
        
        friend PopulationCaller;
        
        Latents(const std::vector<SampleName>& samples,
                const std::vector<Haplotype>&,
                std::vector<Genotype<Haplotype>>&& genotypes,
                ModelInferences&&);
        Latents(const std::vector<SampleName>& samples,
                const std::vector<Haplotype>&,
                std::vector<Genotype<Haplotype>>&& genotypes,
                ModelInferences&&, ModelInferences&&);
        
        std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors() const noexcept override;
        std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors() const noexcept override;
        
    private:
        std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors_;
        std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors_;
        
        boost::optional<ModelInferences> dummy_latents_;
        
        double model_log_evidence_;
        
        HaplotypeProbabilityMap
        calculate_haplotype_posteriors(const std::vector<Haplotype>& haplotypes);
    };
    
    Parameters parameters_;
    
    CallTypeSet do_get_call_types() const override;
    
    std::unique_ptr<Caller::Latents>
    infer_latents(const std::vector<Haplotype>& haplotypes,
                  const HaplotypeLikelihoodCache& haplotype_likelihoods) const override;
    
    std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, const Caller::Latents& latents) const override;
    
    std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, const Latents& latents) const;
    
    std::vector<std::unique_ptr<ReferenceCall>>
    call_reference(const std::vector<Allele>& alleles, const Caller::Latents& latents,
                   const ReadMap& reads) const override;
};
} // namespace octopus

#endif /* defined(__Octopus__population_caller__) */
