//
//  individual_caller.hpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef individual_caller_hpp
#define individual_caller_hpp

#include <vector>
#include <string>
#include <memory>

#include <boost/optional.hpp>

#include "common.hpp"
#include "variant_caller.hpp"
#include "individual_genotype_model.hpp"
#include "variant_call.hpp"

class GenomicRegion;
class ReadPipe;
class Variant;
class HaplotypeLikelihoodCache;

namespace Octopus
{
class IndividualVariantCaller : public VariantCaller
{
public:
    using VariantCaller::CallTypeSet;
    using VariantCaller::CallerComponents;
    
    struct CallerParameters
    {
        double min_variant_posterior;
        double min_refcall_posterior;
        unsigned ploidy;
    };
    
    IndividualVariantCaller() = delete;
    
    IndividualVariantCaller(CallerComponents&& components, VariantCaller::CallerParameters general_parameters,
                            CallerParameters specific_parameters);
    
    ~IndividualVariantCaller() = default;

    IndividualVariantCaller(const IndividualVariantCaller&)            = delete;
    IndividualVariantCaller& operator=(const IndividualVariantCaller&) = delete;
    IndividualVariantCaller(IndividualVariantCaller&&)                 = delete;
    IndividualVariantCaller& operator=(IndividualVariantCaller&&)      = delete;
    
private:
    class Latents;
    
    double min_variant_posterior_;
    double min_refcall_posterior_;
    unsigned ploidy_;
    
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
    
    std::vector<std::unique_ptr<ReferenceCall>>
    call_reference(const std::vector<Allele>& alleles, const Latents& latents,
                   const ReadMap& reads) const;
    
    const SampleIdType& sample() const noexcept;
};

class IndividualVariantCaller::Latents : public CallerLatents
{
public:
    using ModelInferences = GenotypeModel::Individual::InferredLatents;
    
    using CallerLatents::HaplotypeProbabilityMap;
    using CallerLatents::GenotypeProbabilityMap;
    
    friend IndividualVariantCaller;
    
    Latents() = delete;
    
    Latents(const SampleIdType& sample, const std::vector<Haplotype>& haplotypes,
            std::vector<Genotype<Haplotype>>&& genotypes, ModelInferences&& latents);
    
    std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors() const noexcept override;
    std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors() const noexcept override;
    
private:
    std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors_;
    std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors_;
    
    double model_log_evidence_;
    
    HaplotypeProbabilityMap calculate_haplotype_posteriors(const std::vector<Haplotype>& haplotypes);
};
} // namespace Octopus

#endif /* individual_caller_hpp */
