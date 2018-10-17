// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef cell_caller_hpp
#define cell_caller_hpp

#include <vector>
#include <string>
#include <memory>
#include <functional>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "basics/phred.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/models/mutation/coalescent_model.hpp"
#include "core/models/genotype/genotype_prior_model.hpp"
#include "caller.hpp"

namespace octopus {

class GenomicRegion;
class ReadPipe;
class Variant;
class HaplotypeLikelihoodArray;
class VariantCall;

class CellCaller : Caller
{
public:
    using Caller::CallTypeSet;
    
    struct Parameters
    {
        unsigned ploidy;
        boost::optional<CoalescentModel::Parameters> prior_model_params;
        Phred<double> min_variant_posterior, min_refcall_posterior;
        bool deduplicate_haplotypes_with_prior_model = false;
    };
    
    CellCaller() = delete;
    
    CellCaller(Caller::Components&& components,
               Caller::Parameters general_parameters,
               Parameters specific_parameters);
    
    CellCaller(const CellCaller&)            = delete;
    CellCaller& operator=(const CellCaller&) = delete;
    CellCaller(CellCaller&&)                 = delete;
    CellCaller& operator=(CellCaller&&)      = delete;
    
    ~CellCaller() = default;

private:
    class Latents;
    friend Latents;
    
    Parameters parameters_;
    
    std::string do_name() const override;
    CallTypeSet do_call_types() const override;
    unsigned do_min_callable_ploidy() const override;
    unsigned do_max_callable_ploidy() const override;
    
    std::size_t do_remove_duplicates(std::vector<Haplotype>& haplotypes) const override;
    
    std::unique_ptr<Caller::Latents>
    infer_latents(const std::vector<Haplotype>& haplotypes,
                  const HaplotypeLikelihoodArray& haplotype_likelihoods) const override;
    
    boost::optional<double>
    calculate_model_posterior(const std::vector<Haplotype>& haplotypes,
                              const HaplotypeLikelihoodArray& haplotype_likelihoods,
                              const Caller::Latents& latents) const override;
    
    boost::optional<double>
    calculate_model_posterior(const std::vector<Haplotype>& haplotypes,
                              const HaplotypeLikelihoodArray& haplotype_likelihoods,
                              const Latents& latents) const;
    
    std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, const Caller::Latents& latents) const override;
    
    std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, const Latents& latents) const;
    
    std::vector<std::unique_ptr<ReferenceCall>>
    call_reference(const std::vector<Allele>& alleles, const Caller::Latents& latents,
                   const ReadPileupMap& pileup) const override;
    
    std::vector<std::unique_ptr<ReferenceCall>>
    call_reference(const std::vector<Allele>& alleles, const Latents& latents,
                   const ReadPileupMap& pileup) const;
};

class CellCaller::Latents : public Caller::Latents
{
public:
    using Caller::Latents::HaplotypeProbabilityMap;
    using Caller::Latents::GenotypeProbabilityMap;
    
    Latents() = delete;
    
    std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors() const noexcept override;
    std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors() const noexcept override;

private:
    mutable std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors_;
    mutable std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors_;
    
    friend CellCaller;
};

} // namespace octopus

#endif
