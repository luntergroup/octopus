// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef family_caller_hpp
#define family_caller_hpp

#include <vector>
#include <string>

#include "caller.hpp"
#include "basics/pedigree.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/models/mutation/coalescent_model.hpp"
#include "core/models/mutation/denovo_model.hpp"
#include "core/models/genotype/pedigree_prior_model.hpp"
#include "core/models/genotype/family_model.hpp"

namespace octopus {

class GenomicRegion;
class Variant;
class VcfRecord;
class ReadPipe;

class FamilyCaller : public Caller
{
public:
    struct Parameters
    {
        Pedigree family;
        std::vector<unsigned> ploidies;
        boost::optional<CoalescentModel::Parameters> population_prior_model_params;
        DeNovoModel::Parameters denovo_model_params;
        Phred<double> min_variant_posterior, min_denovo_posterior, min_refcall_posterior;
        boost::optional<std::size_t> max_genotype_combinations;
        bool deduplicate_haplotypes_with_germline_model = true;
    };
    
    FamilyCaller() = delete;
    
    FamilyCaller(Caller::Components&& components,
                 Caller::Parameters general_parameters,
                 Parameters specific_parameters);
    
    FamilyCaller(const FamilyCaller&)            = delete;
    FamilyCaller& operator=(const FamilyCaller&) = delete;
    FamilyCaller(FamilyCaller&&)                 = delete;
    FamilyCaller& operator=(FamilyCaller&&)      = delete;
    
    ~FamilyCaller() = default;
    
private:
    class Latents;
    
    using GenotypeBlock = MappableBlock<Genotype<IndexedHaplotype<>>>;
    
    Parameters parameters_;
    std::vector<unsigned> unique_ploidies_;
    std::vector<std::pair<std::size_t, std::pair<std::size_t, std::size_t>>> offspring_with_two_parents_;
    std::vector<std::pair<std::size_t, std::size_t>> offspring_with_one_parent_;
    
    std::string do_name() const override;
    CallTypeSet do_call_types() const override;
    unsigned do_min_callable_ploidy() const override;
    unsigned do_max_callable_ploidy() const override;
    
    std::size_t do_remove_duplicates(HaplotypeBlock& haplotypes) const override;
    
    std::unique_ptr<Caller::Latents>
    infer_latents(const HaplotypeBlock& haplotypes,
                  const HaplotypeLikelihoodArray& haplotype_likelihoods) const override;
    
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

    std::unique_ptr<PedigreePriorModel> make_prior_model(const HaplotypeBlock& haplotypes) const;
};

class FamilyCaller::Latents : public Caller::Latents
{
public:
    using IndexedHaplotypeBlock = MappableBlock<IndexedHaplotype<>>;
    using Genotypeblock = MappableBlock<Genotype<IndexedHaplotype<>>>;
    using ModelInferences = model::FamilyModel::InferredLatents;
    friend FamilyCaller;

    Latents(const std::vector<SampleName>& samples,
            const Pedigree& family,
            IndexedHaplotypeBlock haplotypes,
            Genotypeblock genotypes,
            ModelInferences latents);
    
    std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors() const noexcept override;
    std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors() const noexcept override;

private:
    Genotypeblock genotypes_;
    IndexedHaplotypeBlock haplotypes_;
    ModelInferences inferences_;
    std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors_;
    std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors_;
};

} // namespace octopus

#endif
