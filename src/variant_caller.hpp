//
//  variant_caller.hpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_variant_caller_hpp
#define Octopus_variant_caller_hpp

#include <vector>
#include <string>
#include <functional>
#include <memory>
#include <deque>

#include "common.hpp"
#include "reference_genome.hpp"
#include "read_pipe.hpp"
#include "candidate_variant_generator.hpp"
#include "haplotype_prior_model.hpp"
#include "haplotype_likelihood_cache.hpp"
#include "probability_matrix.hpp"
#include "phaser.hpp"
#include "vcf_record.hpp"
#include "progress_meter.hpp"

class GenomicRegion;
class Variant;
class Haplotype;

namespace Octopus
{
class VariantCaller
{
public:
    enum class RefCallType { None, Blocked, Positional };
    
    struct CallerParameters
    {
        CallerParameters() = default;
        explicit CallerParameters(unsigned max_haplotypes, RefCallType refcall_type,
                                  bool call_sites_only);
        ~CallerParameters() = default;
        
        unsigned max_haplotypes;
        RefCallType refcall_type;
        bool call_sites_only;
    };
    
    using ReadMap = Octopus::ReadMap;
    
    VariantCaller() = delete;
    
    explicit VariantCaller(const ReferenceGenome& reference,
                           ReadPipe& read_pipe,
                           CandidateVariantGenerator&& candidate_generator,
                           std::unique_ptr<HaplotypePriorModel> haplotype_prior_model,
                           CallerParameters parameters);
    
    virtual ~VariantCaller() = default;
    
    VariantCaller(const VariantCaller&)            = delete;
    VariantCaller& operator=(const VariantCaller&) = delete;
    VariantCaller(VariantCaller&&)                 = delete;
    VariantCaller& operator=(VariantCaller&&)      = delete;
    
    std::deque<VcfRecord> call_variants(const GenomicRegion& call_region,
                                        ProgressMeter& progress_meter) const;
    
protected:
    using HaplotypeReference = std::reference_wrapper<const Haplotype>;
    using HaplotypePriorMap  = std::unordered_map<HaplotypeReference, double>;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    std::reference_wrapper<ReadPipe> read_pipe_;
    
    const RefCallType refcall_type_ = RefCallType::Positional;
    const bool call_sites_only_ = false;
    
    bool refcalls_requested() const noexcept;
    
    struct CallerLatents
    {
        using HaplotypePosteriorMap = std::unordered_map<HaplotypeReference, double>;
        using GenotypePosteriorMap  = ProbabilityMatrix<Genotype<Haplotype>>;
        
        // Return shared_ptr as the caller may not actually need to use these latents itself,
        // they just need to be constructable on demand. But if the caller DOES use them, then
        // we avoid copying.
        virtual std::shared_ptr<HaplotypePosteriorMap> get_haplotype_posteriors() const = 0;
        virtual std::shared_ptr<GenotypePosteriorMap> get_genotype_posteriors() const = 0;
        
        virtual ~CallerLatents() = default;
    };
    
private:
    mutable CandidateVariantGenerator candidate_generator_;
    
    unsigned max_haplotypes_;
    
    std::unique_ptr<HaplotypePriorModel> haplotype_prior_model_;
    
    bool done_calling(const GenomicRegion& region) const noexcept;
    
    virtual std::unique_ptr<CallerLatents>
    infer_latents(const std::vector<SampleIdType>& samples,
                  const std::vector<Haplotype>& haplotypes,
                  const HaplotypePriorMap& haplotype_priors,
                  const HaplotypeLikelihoodCache& haplotype_likelihoods) const = 0;
    
    virtual std::vector<VcfRecord::Builder>
    call_variants(const std::vector<Variant>& candidates,
                  const std::vector<Allele>& callable_alleles,
                  CallerLatents* latents,
                  const Phaser::PhaseSet& phase_set,
                  const ReadMap& reads) const = 0;
    
    std::vector<Haplotype>
    get_removable_haplotypes(const std::vector<Haplotype>& haplotypes,
                             const CallerLatents::HaplotypePosteriorMap& haplotype_posteriors,
                             const GenomicRegion& region) const;
};

std::vector<Allele>
generate_callable_alleles(const GenomicRegion& region, const std::vector<Variant>& variants,
                          VariantCaller::RefCallType refcall_type, const ReferenceGenome& reference);

std::vector<Allele>
generate_candidate_reference_alleles(const std::vector<Allele>& callable_alleles,
                                     const std::vector<GenomicRegion>& called_regions,
                                     const std::vector<Variant>& candidates,
                                     VariantCaller::RefCallType refcall_type);
} // namespace Octopus

#endif
