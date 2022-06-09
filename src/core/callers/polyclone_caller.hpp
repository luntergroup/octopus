// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef polyclone_caller_hpp
#define polyclone_caller_hpp

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
#include "core/models/genotype/individual_model.hpp"
#include "core/models/genotype/subclone_model.hpp"
#include "utils/maths.hpp"
#include "caller.hpp"

namespace octopus {

class GenomicRegion;
class ReadPipe;
class Variant;
class HaplotypeLikelihoodArray;
class VariantCall;

class PolycloneCaller : public Caller
{
public:
    using Caller::CallTypeSet;
    
    struct Parameters
    {
        boost::optional<CoalescentModel::Parameters> prior_model_params;
        Phred<double> min_variant_posterior;
        bool deduplicate_haplotypes_with_germline_model = false;
        unsigned max_clones = 3;
        boost::optional<std::size_t> max_genotypes = 10'000;
        boost::optional<unsigned> max_vb_seeds = boost::none; // Use default if none
        std::function<double(unsigned)> clonality_prior = [] (unsigned clonality) { return maths::geometric_pdf(clonality, 0.99); };
        double clone_mixture_prior_concentration = 1;
        bool haplogroup_prior = true;
    };
    
    PolycloneCaller() = delete;
    
    PolycloneCaller(Caller::Components&& components,
                    Caller::Parameters general_parameters,
                    Parameters specific_parameters);
    
    PolycloneCaller(const PolycloneCaller&)            = delete;
    PolycloneCaller& operator=(const PolycloneCaller&) = delete;
    PolycloneCaller(PolycloneCaller&&)                 = delete;
    PolycloneCaller& operator=(PolycloneCaller&&)      = delete;
    
    ~PolycloneCaller() = default;

private:
    class Latents;
    friend Latents;
    
    Parameters parameters_;
    
    struct ModelProbabilities
    {
        double clonal, subclonal;
    };
    
    using IndexedHaplotypeBlock = MappableBlock<IndexedHaplotype<>>;
    using GenotypeBlock = MappableBlock<Genotype<IndexedHaplotype<>>>;
    
    std::string do_name() const override;
    CallTypeSet do_call_types() const override;
    unsigned do_min_callable_ploidy() const override;
    unsigned do_max_callable_ploidy() const override;
    
    std::size_t do_remove_duplicates(HaplotypeBlock& haplotypes) const override;
    
    std::unique_ptr<Caller::Latents>
    infer_latents(const HaplotypeBlock& haplotypes,
                  const HaplotypeLikelihoodArray& haplotype_likelihoods,
                  OptionalThreadPool workers) const override;
    
    std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, const Caller::Latents& latents, OptionalThreadPool workers) const override;
    
    std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, const Latents& latents) const;
    
    std::vector<std::unique_ptr<ReferenceCall>>
    call_reference(const std::vector<Allele>& alleles, const Caller::Latents& latents,
                   const ReadPileupMap& pileup) const override;
    
    std::vector<std::unique_ptr<ReferenceCall>>
    call_reference(const std::vector<Allele>& alleles, const Latents& latents,
                   const ReadPileupMap& pileup) const;
    
    const SampleName& sample() const noexcept;
    void fit_sublone_model(const HaplotypeBlock& haplotypes,
                           const IndexedHaplotypeBlock& indexed_haplotypes,
                           const HaplotypeLikelihoodArray& haplotype_likelihoods,
                           GenotypePriorModel& genotype_prior_model,
                           const model::IndividualModel::InferredLatents& haploid_latents,
                           GenotypeBlock& prev_genotypes,
                           model::SubcloneModel::InferredLatents& sublonal_inferences,
                           OptionalThreadPool workers) const;
    void fit_haplogroup_model(const HaplotypeBlock& haplotypes,
                              const IndexedHaplotypeBlock& indexed_haplotypes,
                              const HaplotypeLikelihoodArray& haplotype_likelihoods,
                              const GenotypeBlock& polyploid_genotypes,
                              const GenotypePriorModel& genotype_prior_model,
                              const model::SubcloneModel::InferredLatents& sublonal_inferences) const;
    
    std::unique_ptr<GenotypePriorModel> make_prior_model(const HaplotypeBlock& haplotypes) const;
    
    // debug
    void log(const Latents& latents) const;
};

class PolycloneCaller::Latents : public Caller::Latents
{
public:
    using GenotypeBlock = MappableBlock<Genotype<IndexedHaplotype<>>>;
    
    using HaploidModelInferences = model::IndividualModel::InferredLatents;
    using SubloneModelInferences = model::SubcloneModel::InferredLatents;
    
    using Caller::Latents::HaplotypeProbabilityMap;
    using Caller::Latents::GenotypeProbabilityMap;
    
    Latents() = delete;
    
    Latents(GenotypeBlock haploid_genotypes,
            GenotypeBlock polyploid_genotypes,
            HaploidModelInferences haploid_model_inferences,
            SubloneModelInferences subclone_model_inferences,
            const SampleName& sample,
            const std::function<double(unsigned)>& clonality_prior);
    
    std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors() const noexcept override;
    std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors() const noexcept override;

private:
    GenotypeBlock haploid_genotypes_, polyploid_genotypes_;
    HaploidModelInferences haploid_model_inferences_;
    SubloneModelInferences subclone_model_inferences_;
    PolycloneCaller::ModelProbabilities model_log_posteriors_;
    SampleName sample_;
    
    mutable std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors_, genotype_log_posteriors_;
    mutable std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors_;
    
    friend PolycloneCaller;
};

} // namespace octopus

#endif
