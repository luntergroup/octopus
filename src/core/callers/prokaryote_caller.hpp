// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef prokaryote_caller_hpp
#define prokaryote_caller_hpp

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
class HaplotypeLikelihoodCache;
class VariantCall;

class ProkaryoteCaller : public Caller
{
public:
    using Caller::CallTypeSet;
    
    struct Parameters
    {
        boost::optional<CoalescentModel::Parameters> prior_model_params;
        Phred<double> min_variant_posterior, min_refcall_posterior;
        bool deduplicate_haplotypes_with_germline_model = false;
        unsigned max_clones = 3;
        std::function<double(unsigned)> clonality_prior = [] (unsigned clonality) { return maths::geometric_pdf(clonality, 0.5); };
    };
    
    ProkaryoteCaller() = delete;
    
    ProkaryoteCaller(Caller::Components&& components,
                     Caller::Parameters general_parameters,
                     Parameters specific_parameters);
    
    ProkaryoteCaller(const ProkaryoteCaller&)            = delete;
    ProkaryoteCaller& operator=(const ProkaryoteCaller&) = delete;
    ProkaryoteCaller(ProkaryoteCaller&&)                 = delete;
    ProkaryoteCaller& operator=(ProkaryoteCaller&&)      = delete;
    
    ~ProkaryoteCaller() = default;

private:
    class Latents;
    friend Latents;
    
    Parameters parameters_;
    
    struct ModelProbabilities
    {
        double clonal, subclonal;
    };
    
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
                   const ReadMap& reads) const override;
    
    std::vector<std::unique_ptr<ReferenceCall>>
    call_reference(const std::vector<Allele>& alleles, const Latents& latents,
                   const ReadMap& reads) const;
    
    const SampleName& sample() const noexcept;
    
    std::unique_ptr<GenotypePriorModel> make_prior_model(const std::vector<Haplotype>& haplotypes) const;
};

class ProkaryoteCaller::Latents : public Caller::Latents
{
public:
    using HaploidModelInferences = model::IndividualModel::InferredLatents;
    using SubloneModelInferences = model::SubcloneModel::InferredLatents;
    
    using Caller::Latents::HaplotypeProbabilityMap;
    using Caller::Latents::GenotypeProbabilityMap;
    
    Latents() = delete;
    
    Latents(std::vector<Genotype<Haplotype>> haploid_genotypes, std::vector<Genotype<Haplotype>> polyploid_genotypes,
            HaploidModelInferences haploid_model_inferences, SubloneModelInferences subclone_model_inferences,
            const SampleName& sample, const std::function<double(unsigned)>& clonality_prior);
    
    std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors() const noexcept override;
    std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors() const noexcept override;

private:
    std::vector<Genotype<Haplotype>> haploid_genotypes_, polyploid_genotypes_;
    HaploidModelInferences haploid_model_inferences_;
    SubloneModelInferences subclone_model_inferences_;
    ProkaryoteCaller::ModelProbabilities model_posteriors_;
    SampleName sample_;
    
    mutable std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors_;
    mutable std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors_;
    
    friend ProkaryoteCaller;
};

} // namespace octopus

#endif
