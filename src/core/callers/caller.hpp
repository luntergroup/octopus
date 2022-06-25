// Copyright (c) 2015-2021 Daniel Cooke
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
#include <boost/variant.hpp>

#include "config/common.hpp"
#include "basics/genomic_region.hpp"
#include "basics/read_pileup.hpp"
#include "core/types/variant.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/indexed_haplotype.hpp"
#include "core/tools/coretools.hpp"
#include "core/models/haplotype_likelihood_array.hpp"
#include "core/tools/vcf_record_factory.hpp"
#include "containers/mappable_flat_set.hpp"
#include "containers/probability_matrix.hpp"
#include "containers/mappable_block.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/reference/reference_genome.hpp"
#include "readpipe/read_pipe.hpp"
#include "utils/memory_footprint.hpp"
#include "logging/progress_meter.hpp"
#include "logging/logging.hpp"

namespace octopus {

class Call;
struct CallWrapper;
class VariantCall;
class ReferenceCall;
class ThreadPool;

class Caller
{
public:
    using CallTypeSet = std::set<std::type_index>;
    
    enum class RefCallType { none, blocked, positional };
    enum class ModelPosteriorPolicy { all, off, special };
    
    struct Components
    {
        std::reference_wrapper<const ReferenceGenome> reference;
        std::reference_wrapper<const ReadPipe> read_pipe;
        VariantGenerator candidate_generator;
        HaplotypeGenerator::Builder haplotype_generator_builder;
        HaplotypeLikelihoodModel likelihood_model;
        Phaser phaser;
        boost::optional<BadRegionDetector> bad_region_detector = boost::none;
    };
    
    struct Parameters
    {
        RefCallType refcall_type;
        boost::optional<Phred<double>> refcall_block_merge_threshold, max_refcall_posterior;
        bool call_sites_only;
        unsigned max_haplotypes;
        double haplotype_extension_threshold, saturation_limit;
        ModelPosteriorPolicy model_posterior_policy;
        bool protect_reference_haplotype;
        boost::optional<MemoryFootprint> target_max_memory;
        ExecutionPolicy execution_policy;
        ReadLinkageType read_linkage;
        bool try_early_phase_detection;
    };
    
    using ReadMap = octopus::ReadMap;

    using OptionalThreadPool = boost::optional<ThreadPool&>;
    
    Caller() = delete;
    
    Caller(Components&& components, Parameters parameters);
    
    Caller(const Caller&)            = delete;
    Caller& operator=(const Caller&) = delete;
    Caller(Caller&&)                 = delete;
    Caller& operator=(Caller&&)      = delete;
    
    virtual ~Caller() = default;
    
    std::string name() const;
    
    CallTypeSet call_types() const;
    
    unsigned min_callable_ploidy() const;
    unsigned max_callable_ploidy() const;
    
    std::deque<VcfRecord> 
    call(const GenomicRegion& call_region, 
        ProgressMeter& progress_meter,
        OptionalThreadPool workers = boost::none) const;
    
    std::vector<VcfRecord> 
    regenotype(const std::vector<Variant>& variants, 
               ProgressMeter& progress_meter) const;
    
protected:
    using HaplotypeBlock = HaplotypeGenerator::HaplotypeBlock;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    
    std::vector<SampleName> samples_;
    
    mutable boost::optional<logging::DebugLogger> debug_log_;
    mutable boost::optional<logging::TraceLogger> trace_log_;
    
    struct Latents
    {
        using HaplotypeProbabilityMap = std::unordered_map<IndexedHaplotype<>, double>;
        using GenotypeProbabilityMap  = ProbabilityMatrix<Genotype<IndexedHaplotype<>>>;
        
        virtual ~Latents() noexcept = default;
        
        // Return shared_ptr as the caller may not actually need to use these latents itself,
        // they just need to be constructable on demand. But if the caller does use them, then
        // we avoid copying.
        virtual std::shared_ptr<HaplotypeProbabilityMap> haplotype_posteriors() const = 0;
        virtual std::shared_ptr<GenotypeProbabilityMap> genotype_posteriors() const = 0;
    };
    
private:
    enum class GeneratorStatus { good, skipped, done };
    
    using GenotypeCallMap = Phaser::GenotypeCallMap;
    
    std::reference_wrapper<const ReadPipe> read_pipe_;
    mutable VariantGenerator candidate_generator_;
    HaplotypeGenerator::Builder haplotype_generator_builder_;
    HaplotypeLikelihoodModel likelihood_model_;
    Phaser phaser_;
    boost::optional<BadRegionDetector> bad_region_detector_;
    Parameters parameters_;
    
    // virtual methods
    
    virtual std::string do_name() const = 0;
    virtual CallTypeSet do_call_types() const = 0;
    virtual unsigned do_min_callable_ploidy() const { return 1; }
    virtual unsigned do_max_callable_ploidy() const { return min_callable_ploidy(); };

protected:
    virtual std::size_t do_remove_duplicates(HaplotypeBlock& haplotypes) const;
    
    struct ModelPosterior
    {
        using Probability = double;
        std::vector<boost::optional<Probability>> samples;
        boost::optional<Probability> joint;
    };

    using ReadPileupMap = std::unordered_map<SampleName, ReadPileups>;
    
    boost::optional<MemoryFootprint> target_max_memory() const noexcept;
    ExecutionPolicy exucution_policy() const noexcept;

private:
    virtual std::unique_ptr<Latents>
    infer_latents(const HaplotypeBlock& haplotypes,
                  const HaplotypeLikelihoodArray& haplotype_likelihoods,
                  OptionalThreadPool workers = boost::none) const = 0;
    
    virtual Genotype<IndexedHaplotype<>> call_genotype(const Latents& latents, const SampleName& sample) const;
    
    virtual boost::optional<ModelPosterior>
    calculate_model_posterior(const HaplotypeBlock& haplotypes,
                              const HaplotypeLikelihoodArray& haplotype_likelihoods,
                              const Latents& latents) const { return boost::none; }
    
    virtual std::vector<std::unique_ptr<VariantCall>>
    call_variants(const std::vector<Variant>& candidates, 
                  const Latents& latents,
                  OptionalThreadPool workers = boost::none) const = 0;
    
    virtual std::vector<std::unique_ptr<ReferenceCall>>
    call_reference(const std::vector<Allele>& alleles, const Latents& latents,
                   const ReadPileupMap& pileups) const = 0;
    
    // helper methods
    
    boost::optional<TemplateMap> make_read_templates(const ReadMap& reads) const;
    std::deque<CallWrapper>
    call_variants(const GenomicRegion& call_region,
                  const MappableFlatSet<Variant>& candidates,
                  const ReadMap& reads,
                  const boost::optional<TemplateMap>& read_templates,
                  HaplotypeGenerator& haplotype_generator,
                  ProgressMeter& progress_meter,
                  OptionalThreadPool workers) const;
    bool refcalls_requested() const noexcept;
    MappableFlatSet<Variant> 
    generate_candidate_variants(const GenomicRegion& region, OptionalThreadPool workers) const;
    HaplotypeGenerator 
    make_haplotype_generator(const MappableFlatSet<Variant>& candidates,
                             const ReadMap& reads,
                             const boost::optional<TemplateMap>& read_templates) const;
    HaplotypeLikelihoodArray make_haplotype_likelihood_cache() const;
    VcfRecordFactory make_record_factory(const ReadMap& reads) const;
    std::vector<Haplotype>
    filter(HaplotypeBlock& haplotypes, const HaplotypeLikelihoodArray& haplotype_likelihoods,
           const std::deque<Haplotype>& protected_haplotypes) const;
    bool compute_haplotype_likelihoods(HaplotypeLikelihoodArray& haplotype_likelihoods, const GenomicRegion& active_region,
                                       const HaplotypeBlock& haplotypes, const MappableFlatSet<Variant>& candidates,
                                       const boost::variant<ReadMap, TemplateMap>& active_reads,
                                       OptionalThreadPool workers) const;
    std::vector<std::reference_wrapper<const Haplotype>>
    get_removable_haplotypes(const HaplotypeBlock& haplotypes, const HaplotypeLikelihoodArray& haplotype_likelihoods,
                             const Latents::HaplotypeProbabilityMap& haplotype_posteriors,
                             const std::deque<Haplotype>& protected_haplotypes, unsigned max_to_remove) const;
    GeneratorStatus
    generate_active_haplotypes(const GenomicRegion& call_region, HaplotypeGenerator& haplotype_generator,
                               GenomicRegion& active_region, boost::optional<GenomicRegion>& next_active_region,
                               HaplotypeBlock& haplotypes, HaplotypeBlock& next_haplotypes, 
                               boost::optional<GenomicRegion> backtrack_region) const;
    GeneratorStatus
    generate_next_active_haplotypes(HaplotypeBlock& next_haplotypes,
                                    boost::optional<GenomicRegion>& next_active_region,
                                    boost::optional<GenomicRegion>& backtrack_region,
                                    HaplotypeGenerator& haplotype_generator) const;
    void remove_duplicates(HaplotypeBlock& haplotypes) const;
    bool filter_haplotypes(HaplotypeBlock& haplotypes, HaplotypeGenerator& haplotype_generator,
                           HaplotypeLikelihoodArray& haplotype_likelihoods,
                           const std::deque<Haplotype>& protected_haplotypes) const;
    bool is_saturated(const HaplotypeBlock& haplotypes, const Latents& latents) const;
    unsigned count_probable_haplotypes(const Caller::Latents::HaplotypeProbabilityMap& haplotype_posteriors) const;
    void filter_haplotypes(bool prefilter_had_removal_impact, const HaplotypeBlock& haplotypes,
                           HaplotypeGenerator& haplotype_generator, const HaplotypeLikelihoodArray& haplotype_likelihoods,
                           const Latents& latents, const std::deque<Haplotype>& protected_haplotypes) const;
    bool try_early_detect_phase_regions(const MappableBlock<Haplotype>& haplotypes,
                                        const MappableFlatSet<Variant>& candidates,
                                        const GenomicRegion& active_region,
                                        const Latents& latents,
                                        const boost::optional<GenomicRegion>& backtrack_region) const;
    boost::optional<GenomicRegion>
    find_phased_head(const MappableBlock<Haplotype>& haplotypes, const MappableFlatSet<Variant>& candidates,
                     const GenomicRegion& active_region, const Latents& latents) const;
    void call_variants(const GenomicRegion& active_region, const GenomicRegion& call_region,
                       const boost::optional<GenomicRegion>& next_active_region,
                       const boost::optional<GenomicRegion>& backtrack_region,
                       const MappableFlatSet<Variant>& candidates, const HaplotypeBlock& haplotypes,
                       const HaplotypeLikelihoodArray& haplotype_likelihoods, const ReadMap& reads,
                       const Latents& latents, std::deque<CallWrapper>& result,
                       boost::optional<GenomicRegion>& prev_called_region, GenomicRegion& completed_region) const;
    GenotypeCallMap get_genotype_calls(const Latents& latents) const;
    std::deque<Haplotype> get_called_haplotypes(const Latents& latents) const;
    void set_model_posteriors(std::vector<CallWrapper>& calls, const Latents& latents,
                              const HaplotypeBlock& haplotypes,
                              const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    void set_phasing(std::vector<CallWrapper>& calls, const Latents& latents,
                     const HaplotypeBlock& haplotypes, const GenomicRegion& call_region) const;
    bool done_calling(const GenomicRegion& region) const noexcept;
    bool is_merge_block_refcalling() const noexcept;
    std::vector<CallWrapper> call_reference(const GenomicRegion& region, const ReadMap& reads) const;
    std::vector<CallWrapper> call_reference_helper(const std::vector<Allele>& alleles, const Latents& latents,
                                                   const ReadPileupMap& pileups) const;
    std::vector<Allele>
    generate_reference_alleles(const GenomicRegion& region,
                               const std::vector<CallWrapper>& calls) const;
    std::vector<Allele> generate_reference_alleles(const GenomicRegion& region) const;
    ReadPileupMap make_pileups(const ReadMap& reads, const Latents& latents, const GenomicRegion& region) const;
    std::vector<std::unique_ptr<ReferenceCall>>
    squash_reference_calls(std::vector<std::unique_ptr<ReferenceCall>> refcalls) const;
};

} // namespace octopus

#endif
