// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef caller_hpp
#define caller_hpp

#include <vector>
#include <string>
#include <functional>
#include <memory>
#include <deque>
#include <typeindex>
#include <set>

#include <boost/optional.hpp>

#include <config/common.hpp>
#include <basics/genomic_region.hpp>
#include <io/reference/reference_genome.hpp>
#include <readpipe/read_pipe.hpp>
#include <core/types/variant.hpp>
#include <core/types/haplotype.hpp>
#include <core/tools/coretools.hpp>
#include <core/models/haplotype_likelihood_cache.hpp>
#include <containers/mappable_flat_set.hpp>
#include <containers/probability_matrix.hpp>
#include <logging/progress_meter.hpp>
#include <logging/logging.hpp>
#include <io/variant/vcf_record.hpp>

#include "utils/vcf_record_factory.hpp"

namespace octopus {

class Call;
class VariantCall;
class ReferenceCall;

class Caller
{
public:
    using CallTypeSet = std::set<std::type_index>;
    
    enum class RefCallType { None, Blocked, Positional };
    
    struct Components;
    struct Parameters;
    
    using ReadMap = octopus::ReadMap;
    
    Caller() = delete;
    
    Caller(Components&& components, Parameters parameters);
    
    Caller(const Caller&)            = delete;
    Caller& operator=(const Caller&) = delete;
    Caller(Caller&&)                 = delete;
    Caller& operator=(Caller&&)      = delete;
    
    virtual ~Caller() = default;
    
    CallTypeSet get_call_types() const;
    
    std::deque<VcfRecord> call(const GenomicRegion& call_region, ProgressMeter& progress_meter) const;
    
    std::vector<VcfRecord> regenotype(const std::vector<Variant>& variants, ProgressMeter& progress_meter) const;
    
protected:
    using HaplotypeReference = std::reference_wrapper<const Haplotype>;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    
    std::vector<SampleName> samples_;
    
    mutable boost::optional<logging::DebugLogger> debug_log_;
    mutable boost::optional<logging::TraceLogger> trace_log_;
    
    struct Latents
    {
        using HaplotypeProbabilityMap = std::unordered_map<HaplotypeReference, double>;
        using GenotypeProbabilityMap  = ProbabilityMatrix<Genotype<Haplotype>>;
        
        virtual ~Latents() noexcept = default;
        
        // Return shared_ptr as the caller may not actually need to use these latents itself,
        // they just need to be constructable on demand. But if the caller does use them, then
        // we avoid copying.
        virtual std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors() const = 0;
        virtual std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors() const = 0;
    };
    
public:
    struct Components
    {
        std::reference_wrapper<const ReferenceGenome> reference;
        std::reference_wrapper<const ReadPipe> read_pipe;
        VariantGenerator candidate_generator;
        HaplotypeGenerator::Builder haplotype_generator_builder;
        Phaser phaser;
    };
    
    struct Parameters
    {
        RefCallType refcall_type;
        bool call_sites_only;
        unsigned max_haplotypes;
        double min_haplotype_posterior;
        bool allow_inactive_flank_scoring;
        bool allow_model_filtering;
    };
    
private:
    std::reference_wrapper<const ReadPipe> read_pipe_;
    
    mutable VariantGenerator candidate_generator_;
    
    HaplotypeGenerator::Builder haplotype_generator_builder_;
    
    Phaser phaser_;
    
    Parameters parameters_;
    
    // virtual methods
    
    virtual CallTypeSet do_get_call_types() const = 0;
    
    virtual std::unique_ptr<Latents>
    infer_latents(const std::vector<Haplotype>& haplotypes,
                  const HaplotypeLikelihoodCache& haplotype_likelihoods) const = 0;
    
    //virtual Genotype<Haplotype> call_genotype(const Latents& latents) const = 0;
    
    virtual boost::optional<double>
    calculate_model_posterior(const std::vector<Haplotype>& haplotypes,
                              const HaplotypeLikelihoodCache& haplotype_likelihoods,
                              const Latents& latents) const { return boost::none; }
    
    virtual std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, const Latents& latents) const = 0;
    
    virtual std::vector<std::unique_ptr<ReferenceCall>>
    call_reference(const std::vector<Allele>& alleles, const Latents& latents,
                   const ReadMap& reads) const = 0;
    
    // helper methods
    
    bool refcalls_requested() const noexcept;
    
    MappableFlatSet<Variant> generate_candidate_variants(const GenomicRegion& region) const;
    
    HaplotypeGenerator make_haplotype_generator(const MappableFlatSet<Variant>& candidates,
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
                             const HaplotypeLikelihoodCache& haplotype_likelihoods,
                             const Latents::HaplotypeProbabilityMap& haplotype_posteriors,
                             unsigned max_to_remove) const;
    
    bool done_calling(const GenomicRegion& region) const noexcept;
    
    std::vector<Allele>
    generate_callable_alleles(const GenomicRegion& region, const std::vector<Variant>& candidates) const;
    
    std::vector<Allele>
    generate_candidate_reference_alleles(const GenomicRegion& region,
                                         const std::vector<Variant>& candidates,
                                         const std::vector<GenomicRegion>& called_regions) const;
};

} // namespace octopus

#endif
