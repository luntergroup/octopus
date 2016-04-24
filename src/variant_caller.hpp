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

#include <boost/optional.hpp>

#include "common.hpp"
#include "reference_genome.hpp"
#include "read_pipe.hpp"
#include "candidate_variant_generator.hpp"
#include "haplotype_likelihood_cache.hpp"
#include "probability_matrix.hpp"
#include "phaser.hpp"
#include "vcf_record.hpp"
#include "call.hpp"
#include "variant_call.hpp"
#include "reference_call.hpp"
#include "progress_meter.hpp"
#include "logging.hpp"

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
                                  bool call_sites_only, bool allow_lagging, double min_phase_score);
        ~CallerParameters() = default;
        
        unsigned max_haplotypes;
        RefCallType refcall_type;
        bool call_sites_only;
        bool lag_haplotype_generation;
        double min_phase_score;
    };
    
    using ReadMap = Octopus::ReadMap;
    
    VariantCaller() = delete;
    
    explicit VariantCaller(const ReferenceGenome& reference,
                           ReadPipe& read_pipe,
                           CandidateVariantGenerator&& candidate_generator,
                           CallerParameters parameters);
    
    virtual ~VariantCaller() = default;
    
    VariantCaller(const VariantCaller&)            = delete;
    VariantCaller& operator=(const VariantCaller&) = delete;
    VariantCaller(VariantCaller&&)                 = delete;
    VariantCaller& operator=(VariantCaller&&)      = delete;
    
    std::deque<VcfRecord> call(const GenomicRegion& call_region, ProgressMeter& progress_meter) const;
    
protected:
    using HaplotypeReference = std::reference_wrapper<const Haplotype>;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    std::reference_wrapper<ReadPipe> read_pipe_;
    
    std::vector<SampleIdType> samples_;
    
    mutable boost::optional<Logging::DebugLogger> debug_log_;
    
    struct CallerLatents
    {
        using HaplotypeProbabilityMap = std::unordered_map<HaplotypeReference, double>;
        using GenotypeProbabilityMap  = ProbabilityMatrix<Genotype<Haplotype>>;
        
        // Return shared_ptr as the caller may not actually need to use these latents itself,
        // they just need to be constructable on demand. But if the caller DOES use them, then
        // we avoid copying.
        virtual std::shared_ptr<HaplotypeProbabilityMap> get_haplotype_posteriors() const = 0;
        virtual std::shared_ptr<GenotypeProbabilityMap> get_genotype_posteriors() const = 0;
        
        virtual ~CallerLatents() = default;
    };
    
private:
    mutable CandidateVariantGenerator candidate_generator_;
    
    unsigned max_haplotypes_;
    double min_haplotype_posterior_;
    bool lag_haplotype_generation_;
    double min_phase_score_;
    
    const RefCallType refcall_type_;
    const bool call_sites_only_;
    
    // virtual methods
    
    virtual std::unique_ptr<CallerLatents>
    infer_latents(const std::vector<Haplotype>& haplotypes,
                  const HaplotypeLikelihoodCache& haplotype_likelihoods) const = 0;
    
    virtual std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, CallerLatents& latents) const = 0;
    
    virtual std::vector<std::unique_ptr<ReferenceCall>>
    call_reference(const std::vector<Allele>& alleles, CallerLatents& latents,
                   const ReadMap& reads) const = 0;
    
    // other private methods
    
    bool done_calling(const GenomicRegion& region) const noexcept;
    bool refcalls_requested() const noexcept;
    
    std::vector<std::reference_wrapper<const Haplotype>>
    get_removable_haplotypes(const std::vector<Haplotype>& haplotypes,
                             const CallerLatents::HaplotypeProbabilityMap& haplotype_posteriors,
                             const GenomicRegion& region) const;
    
    std::vector<Allele>
    generate_callable_alleles(const GenomicRegion& region, const std::vector<Variant>& candidates) const;
    
    std::vector<Allele>
    generate_candidate_reference_alleles(const GenomicRegion& region,
                                         const std::vector<Variant>& candidates,
                                         const std::vector<GenomicRegion>& called_regions) const;
};

} // namespace Octopus

#endif
