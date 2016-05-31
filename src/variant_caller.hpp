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
#include "genomic_region.hpp"
#include "variant.hpp"
#include "haplotype.hpp"
#include "haplotype_generator.hpp"
#include "phaser.hpp"
#include "candidate_variant_generator.hpp"
#include "haplotype_likelihood_cache.hpp"
#include "mappable_flat_set.hpp"
#include "probability_matrix.hpp"
#include "vcf_record.hpp"
#include "vcf_record_factory.hpp"
#include "call.hpp"
#include "variant_call.hpp"
#include "reference_call.hpp"
#include "progress_meter.hpp"
#include "logging.hpp"

namespace Octopus
{
class VariantCaller
{
public:
    enum class RefCallType { None, Blocked, Positional };
    
    struct CallerComponents;
    struct CallerParameters;
    
    using ReadMap = Octopus::ReadMap;
    
    VariantCaller() = delete;
    
    VariantCaller(CallerComponents&& components, CallerParameters parameters);
    
    virtual ~VariantCaller() = default;
    
    VariantCaller(const VariantCaller&)            = delete;
    VariantCaller& operator=(const VariantCaller&) = delete;
    VariantCaller(VariantCaller&&)                 = delete;
    VariantCaller& operator=(VariantCaller&&)      = delete;
    
    std::deque<VcfRecord>
    call(const GenomicRegion& call_region, ProgressMeter& progress_meter) const;
    
//    std::vector<VcfRecord>
//    regenotype(const std::vector<Variant>& variants, ProgressMeter& progress_meter) const;
    
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
        // they just need to be constructable on demand. But if the caller does use them, then
        // we avoid copying.
        virtual std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors() const = 0;
        virtual std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors() const = 0;
        
        virtual ~CallerLatents() noexcept = default;
    };
    
public:
    struct CallerComponents
    {
        CallerComponents() = delete;
        
        CallerComponents(const ReferenceGenome& reference,
                         ReadPipe& read_pipe,
                         CandidateVariantGenerator&& candidate_generator,
                         HaplotypeGenerator::Builder haplotype_generator_builder,
                         Phaser phase);
        
        CallerComponents(const CallerComponents&)            = delete;
        CallerComponents& operator=(const CallerComponents&) = delete;
        CallerComponents(CallerComponents&&)                 = default;
        CallerComponents& operator=(CallerComponents&&)      = default;
        
        std::reference_wrapper<const ReferenceGenome> reference;
        std::reference_wrapper<ReadPipe> read_pipe;
        CandidateVariantGenerator candidate_generator;
        HaplotypeGenerator::Builder haplotype_generator_builder;
        Phaser phaser;
    };
    
    struct CallerParameters
    {
        RefCallType refcall_type;
        bool call_sites_only;
        unsigned max_haplotypes;
        double min_haplotype_posterior;
        bool allow_inactive_flank_scoring;
        bool allow_model_filtering;
    };
    
private:
    mutable CandidateVariantGenerator candidate_generator_;
    
    HaplotypeGenerator::Builder haplotype_generator_builder_;
    
    Phaser phaser_;
    
    CallerParameters parameters_;
    
    // virtual methods
    
    virtual std::unique_ptr<CallerLatents>
    infer_latents(const std::vector<Haplotype>& haplotypes,
                  const HaplotypeLikelihoodCache& haplotype_likelihoods) const = 0;
    
    virtual double
    calculate_dummy_model_posterior(const std::vector<Haplotype>& haplotypes,
                                    const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                    const CallerLatents& latents) const { return 0; }
    
    virtual std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, const CallerLatents& latents) const = 0;
    
    virtual std::vector<std::unique_ptr<ReferenceCall>>
    call_reference(const std::vector<Allele>& alleles, const CallerLatents& latents,
                   const ReadMap& reads) const = 0;
    
    // other private methods
    
    bool refcalls_requested() const noexcept;
    
    MappableFlatSet<Variant> generate_candidates(const GenomicRegion& region) const;
    
    HaplotypeGenerator make_haplotype_generator(const GenomicRegion& region,
                                                const MappableFlatSet<Variant>& candidates,
                                                const ReadMap& reads) const;
    HaplotypeLikelihoodCache make_haplotype_likelihood_cache() const;
    VcfRecordFactory make_record_factory(const ReadMap& reads) const;
    
    std::vector<Haplotype> filter(std::vector<Haplotype>& haplotypes,
                                  const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    
    void populate(HaplotypeLikelihoodCache& haplotype_likelihoods,
                  const GenomicRegion& active_region,
                  const std::vector<Haplotype>& haplotypes, const MappableFlatSet<Variant>& candidates,
                  const ReadMap& active_reads) const;
    
    std::vector<std::reference_wrapper<const Haplotype>>
    get_removable_haplotypes(const std::vector<Haplotype>& haplotypes,
                             const CallerLatents::HaplotypeProbabilityMap& haplotype_posteriors) const;
    
    bool done_calling(const GenomicRegion& region) const noexcept;
    
    std::vector<Allele>
    generate_callable_alleles(const GenomicRegion& region, const std::vector<Variant>& candidates) const;
    
    std::vector<Allele>
    generate_candidate_reference_alleles(const GenomicRegion& region,
                                         const std::vector<Variant>& candidates,
                                         const std::vector<GenomicRegion>& called_regions) const;
};

} // namespace Octopus

#endif
