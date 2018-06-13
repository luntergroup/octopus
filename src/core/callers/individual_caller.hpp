// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef individual_caller_hpp
#define individual_caller_hpp

#include <vector>
#include <string>
#include <memory>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "basics/phred.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/models/mutation/coalescent_model.hpp"
#include "core/models/genotype/genotype_prior_model.hpp"
#include "core/models/genotype/individual_model.hpp"
#include "caller.hpp"

namespace octopus {

class GenomicRegion;
class ReadPipe;
class Variant;
class HaplotypeLikelihoodCache;
class VariantCall;

class IndividualCaller : public Caller
{
public:
    using Caller::CallTypeSet;
    
    struct Parameters
    {
        unsigned ploidy;
        boost::optional<CoalescentModel::Parameters> prior_model_params;
        Phred<double> min_variant_posterior, min_refcall_posterior;
        bool deduplicate_haplotypes_with_germline_model = false;
    };
    
    IndividualCaller() = delete;
    
    IndividualCaller(Caller::Components&& components,
                     Caller::Parameters general_parameters,
                     Parameters specific_parameters);

    IndividualCaller(const IndividualCaller&)            = delete;
    IndividualCaller& operator=(const IndividualCaller&) = delete;
    IndividualCaller(IndividualCaller&&)                 = delete;
    IndividualCaller& operator=(IndividualCaller&&)      = delete;
    
    ~IndividualCaller() = default;
    
private:
    class Latents;
    
    Parameters parameters_;
    
    std::string do_name() const override;
    CallTypeSet do_call_types() const override;
    
    std::size_t do_remove_duplicates(std::vector<Haplotype>& haplotypes) const override;
    
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
                   const ReadPileupMap& pileups) const override;
    
    std::vector<std::unique_ptr<ReferenceCall>>
    call_reference(const std::vector<Allele>& alleles, const Latents& latents,
                   const ReadPileupMap& pileups) const;
    
    const SampleName& sample() const noexcept;
    
    std::unique_ptr<GenotypePriorModel> make_prior_model(const std::vector<Haplotype>& haplotypes) const;
};

class IndividualCaller::Latents : public Caller::Latents
{
public:
    using ModelInferences = model::IndividualModel::InferredLatents;
    
    using Caller::Latents::HaplotypeProbabilityMap;
    using Caller::Latents::GenotypeProbabilityMap;
    
    friend IndividualCaller;
    
    Latents() = delete;
    
    Latents(const SampleName& sample, const std::vector<Haplotype>& haplotypes,
            std::vector<Genotype<Haplotype>>&& genotypes, ModelInferences&& latents);
    
    std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors() const noexcept override;
    std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors() const noexcept override;
    
private:
    std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors_;
    std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors_;
    double model_log_evidence_;
    
    HaplotypeProbabilityMap calculate_haplotype_posteriors(const std::vector<Haplotype>& haplotypes);
};

} // namespace octopus

#endif
