// Copyright (c) 2015-2021 Daniel Cooke
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
class HaplotypeLikelihoodArray;
class VariantCall;

class IndividualCaller : public Caller
{
public:
    using Caller::CallTypeSet;
    
    struct Parameters
    {
        unsigned ploidy;
        boost::optional<CoalescentModel::Parameters> prior_model_params;
        Phred<double> min_variant_posterior;
        bool deduplicate_haplotypes_with_germline_model = false;
        boost::optional<std::size_t> max_genotypes = boost::none;
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
    
    using IndexedHaplotypeBlock = MappableBlock<IndexedHaplotype<>>;
    using GenotypeBlock = MappableBlock<Genotype<IndexedHaplotype<>>>;
    
    Parameters parameters_;
    
    std::string do_name() const override;
    CallTypeSet do_call_types() const override;
    unsigned do_min_callable_ploidy() const override;
    
    std::size_t do_remove_duplicates(HaplotypeBlock& haplotypes) const override;
    
    std::unique_ptr<Caller::Latents>
    infer_latents(const HaplotypeBlock& haplotypes,
                  const HaplotypeLikelihoodArray& haplotype_likelihoods,
                  OptionalThreadPool workers) const override;
    
    boost::optional<ModelPosterior>
    calculate_model_posterior(const HaplotypeBlock& haplotypes,
                              const HaplotypeLikelihoodArray& haplotype_likelihoods,
                              const Caller::Latents& latents) const override;
    
    boost::optional<ModelPosterior>
    calculate_model_posterior(const HaplotypeBlock& haplotypes,
                              const HaplotypeLikelihoodArray& haplotype_likelihoods,
                              const Latents& latents) const;
    
    std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, 
                  const Caller::Latents& latents,
                  OptionalThreadPool workers) const override;
    
    std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, const Latents& latents) const;
    
    std::vector<std::unique_ptr<ReferenceCall>>
    call_reference(const std::vector<Allele>& alleles, const Caller::Latents& latents,
                   const ReadPileupMap& pileups) const override;
    
    std::vector<std::unique_ptr<ReferenceCall>>
    call_reference(const std::vector<Allele>& alleles, const Latents& latents,
                   const ReadPileupMap& pileups) const;
    
    const SampleName& sample() const noexcept;
    
    std::unique_ptr<GenotypePriorModel> make_prior_model(const HaplotypeBlock& haplotypes) const;
    GenotypeBlock 
    propose_genotypes(const HaplotypeBlock& haplotypes,
                      const IndexedHaplotypeBlock& indexed_haplotypes,
                      const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    GenotypeBlock 
    propose_model_check_genotypes(const HaplotypeBlock& haplotypes,
                                  const IndexedHaplotypeBlock& indexed_haplotypes,
                                  const Latents& latents) const;
};

class IndividualCaller::Latents : public Caller::Latents
{
public:
    using IndexedHaplotypeBlock = MappableBlock<IndexedHaplotype<>>;
    
    using ModelInferences = model::IndividualModel::InferredLatents;
    
    using Caller::Latents::HaplotypeProbabilityMap;
    using Caller::Latents::GenotypeProbabilityMap;
    
    friend IndividualCaller;
    
    Latents() = delete;
    
    Latents(const SampleName& sample,
            const IndexedHaplotypeBlock& haplotypes,
            IndividualCaller::GenotypeBlock genotypes,
            ModelInferences&& latents);
    
    std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors() const override;
    std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors() const override;
    std::shared_ptr<GenotypeProbabilityMap> genotype_log_posteriors() const;
    
private:
    mutable std::shared_ptr<GenotypeProbabilityMap> genotype_log_posteriors_, genotype_posteriors_;
    mutable std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors_;
    IndexedHaplotypeBlock haplotypes_;
    IndividualCaller::GenotypeBlock genotypes_;
    ModelInferences latents_;
    SampleName samples_;
    
    HaplotypeProbabilityMap calculate_haplotype_posteriors(const IndexedHaplotypeBlock& haplotypes) const;
};

} // namespace octopus

#endif
