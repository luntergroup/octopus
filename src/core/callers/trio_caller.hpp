// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef trio_caller_hpp
#define trio_caller_hpp

#include <vector>
#include <string>

#include "caller.hpp"
#include "basics/trio.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/models/mutation/coalescent_model.hpp"
#include "core/models/genotype/population_prior_model.hpp"
#include "core/models/genotype/genotype_prior_model.hpp"
#include "core/models/mutation/denovo_model.hpp"
#include "core/models/genotype/trio_model.hpp"

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
        unsigned maternal_ploidy, paternal_ploidy, child_ploidy;
        boost::optional<CoalescentModel::Parameters> germline_prior_model_params;
        DeNovoModel::Parameters denovo_model_params;
        Phred<double> min_variant_posterior, min_denovo_posterior, min_refcall_posterior;
        boost::optional<std::size_t> max_genotype_combinations;
        bool deduplicate_haplotypes_with_germline_model = true;
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
    
    std::string do_name() const override;
    CallTypeSet do_call_types() const override;
    unsigned do_min_callable_ploidy() const override;
    unsigned do_max_callable_ploidy() const override;
    
    std::size_t do_remove_duplicates(HaplotypeBlock& haplotypes) const override;
    
    std::unique_ptr<Caller::Latents>
    infer_latents(const HaplotypeBlock& haplotypes,
                  const HaplotypeLikelihoodArray& haplotype_likelihoods) const override;
    
    boost::optional<double>
    calculate_model_posterior(const HaplotypeBlock& haplotypes,
                              const HaplotypeLikelihoodArray& haplotype_likelihoods,
                              const Caller::Latents& latents) const override;
    
    boost::optional<double>
    calculate_model_posterior(const HaplotypeBlock& haplotypes,
                              const HaplotypeLikelihoodArray& haplotype_likelihoods,
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
    
    std::unique_ptr<PopulationPriorModel> make_prior_model(const HaplotypeBlock& haplotypes) const;
    std::unique_ptr<GenotypePriorModel> make_single_sample_prior_model(const HaplotypeBlock& haplotypes) const;
};

class TrioCaller::Latents : public Caller::Latents
{
public:
    using ModelInferences = model::TrioModel::InferredLatents;
    friend TrioCaller;
    
    Latents(const HaplotypeBlock& haplotypes,
            std::vector<Genotype<Haplotype>>&& genotypes,
            ModelInferences&& latents,
            const Trio& trio);
    Latents(const HaplotypeBlock& haplotypes,
            std::vector<Genotype<Haplotype>>&& maternal_genotypes,
            std::vector<Genotype<Haplotype>>&& paternal_genotypes,
            unsigned child_ploidy,
            ModelInferences&& latents,
            const Trio& trio);
    
    std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors() const noexcept override;
    std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors() const noexcept override;
    
private:
    Trio trio;
    std::vector<Genotype<Haplotype>> maternal_genotypes;
    boost::optional<std::vector<Genotype<Haplotype>>> paternal_genotypes;
    ModelInferences model_latents;
    std::vector<double> marginal_maternal_posteriors, marginal_paternal_posteriors, marginal_child_posteriors;
    mutable std::shared_ptr<GenotypeProbabilityMap> marginal_genotype_posteriors;
    std::shared_ptr<HaplotypeProbabilityMap> marginal_haplotype_posteriors;
    std::vector<Genotype<Haplotype>> concatenated_genotypes_;
    std::vector<double> padded_marginal_maternal_posteriors_, padded_marginal_paternal_posteriors_, padded_marginal_child_posteriors_;
    unsigned child_ploidy_;
    
    void set_genotype_posteriors(const Trio& trio);
    void set_genotype_posteriors_shared_genotypes(const Trio& trio);
    void set_genotype_posteriors_unique_genotypes(const Trio& trio);
    void set_haplotype_posteriors(const HaplotypeBlock& haplotypes);
    void set_haplotype_posteriors_shared_genotypes(const HaplotypeBlock& haplotypes);
    void set_haplotype_posteriors_unique_genotypes(const HaplotypeBlock& haplotypes);
};

} // namespace octopus

#endif
