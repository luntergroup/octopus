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
#include "haplotype_prior_model.hpp"
#include "population_genotype_model.hpp"

class GenomicRegion;
class ReadPipe;
class Variant;
class VcfRecord;
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
                                     std::unique_ptr<HaplotypePriorModel> haplotype_prior_model,
                                     VariantCaller::CallerParameters general_parameters,
                                     CallerParameters specific_parameters);
    
    ~PopulationVariantCaller() = default;
    
    PopulationVariantCaller(const PopulationVariantCaller&)            = delete;
    PopulationVariantCaller& operator=(const PopulationVariantCaller&) = delete;
    PopulationVariantCaller(PopulationVariantCaller&&)                 = delete;
    PopulationVariantCaller& operator=(PopulationVariantCaller&&)      = delete;
    
private:
    using VariantCaller::HaplotypePriorMap;
    
    struct Latents : public CallerLatents
    {
        using ModelLatents = GenotypeModel::Population::Latents;
        
        using CallerLatents::HaplotypePosteiorMap;
        using CallerLatents::GenotypePosteriorMap;
        
        Latents(ModelLatents&&);
        Latents(HaplotypePosteiorMap&&, GenotypePosteriorMap&&);
        
        std::shared_ptr<HaplotypePosteiorMap> get_haplotype_posteriors() const noexcept override
        {
            return haplotype_frequencies_;
        }
        
        std::shared_ptr<GenotypePosteriorMap> get_genotype_posteriors() const noexcept override
        {
            return genotype_posteriors_;
        }
        
        std::shared_ptr<ModelLatents::HaplotypeFrequencyMap> haplotype_frequencies_;
        std::shared_ptr<ModelLatents::GenotypeProbabilityMap> genotype_posteriors_;
    };
    
    mutable GenotypeModel::Population genotype_model_;
    
    const unsigned ploidy_;
    const double min_variant_posterior_ = 0.95;
    const double min_refcall_posterior_ = 0.5;
    
    std::unique_ptr<CallerLatents>
    infer_latents(const std::vector<SampleIdType>& samples,
                  const std::vector<Haplotype>& haplotypes,
                  const HaplotypePriorMap& haplotype_priors,
                  const HaplotypeLikelihoodCache& haplotype_likelihoods) const override;
    
    std::vector<VcfRecord::Builder>
    call_variants(const std::vector<Variant>& candidates,
                  const std::vector<Allele>& callable_alleles,
                  CallerLatents* latents,
                  const Phaser::PhaseSet& phase_set,
                  const ReadMap& reads) const override;
};
} // namespace Octopus

#endif /* defined(__Octopus__population_caller__) */
