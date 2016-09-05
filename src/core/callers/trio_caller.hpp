// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef trio_caller_hpp
#define trio_caller_hpp

#include <vector>
#include <string>

#include "caller.hpp"
#include "core/types/trio.hpp"
#include "core/models/mutation/coalescent_model.hpp"
#include "core/models/mutation/denovo_model.hpp"

namespace octopus {

class GenomicRegion;
class Variant;
class VcfRecord;
class ReadPipe;

class TrioCaller : public Caller
{
public:
    struct Parameters
    {
        Trio trio;
        Phred<double> min_variant_posterior;
        Phred<double> min_refcall_posterior;
        unsigned ploidy;
        CoalescentModel::Parameters germline_prior_model_params;
        DeNovoModel::Parameters denovo_model_params;
    };
    
    TrioCaller() = delete;
    
    TrioCaller(Caller::Components&& components,
               Caller::Parameters general_parameters,
               Parameters specific_parameters);
    
    TrioCaller(const TrioCaller&)            = delete;
    TrioCaller& operator=(const TrioCaller&) = delete;
    TrioCaller(TrioCaller&&)                 = delete;
    TrioCaller& operator=(TrioCaller&&)      = delete;
    
    ~TrioCaller() = default;
    
private:
    class Latents;
    
    Parameters parameters_;
    
    CallTypeSet do_get_call_types() const override;
    
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
    
    std::vector<std::unique_ptr<ReferenceCall>>
    call_reference(const std::vector<Allele>& alleles, const Latents& latents,
                   const ReadMap& reads) const;
};

class TrioCaller::Latents : public Caller::Latents
{
public:
    std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors() const noexcept override { return nullptr; }
    std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors() const noexcept override { return nullptr; }
};
    
} // namespace octopus

#endif
