// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef haplotype_generator_hpp
#define haplotype_generator_hpp

#include <vector>
#include <tuple>
#include <functional>
#include <utility>
#include <stack>
#include <set>
#include <stdexcept>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "basics/genomic_region.hpp"
#include "containers/mappable_flat_set.hpp"
#include "containers/mappable_block.hpp"
#include "core/types/allele.hpp"
#include "readpipe/read_pipe.hpp"
#include "logging/logging.hpp"
#include "genome_walker.hpp"
#include "haplotype_tree.hpp"
#include "core/tools/bad_region_detector.hpp"

namespace octopus {

class Variant;
class Haplotype;

namespace coretools {

/**
    HaplotypeGenerator takes a set of Variant's and turns them into Haplotype's in a controlled
    manner. Haplotype generation always proceeds sequentially (left to right).
 */
class HaplotypeGenerator
{
public:
    struct Policies
    {
        enum class Lagging { none, conservative, moderate, optimistic, aggressive } lagging = Lagging::moderate;
        enum class Extension { minimal, conservative, moderate, aggressive, unlimited } extension = Extension::moderate;
        enum class Backtrack { none, moderate, aggressive } backtrack = Backtrack::moderate;
        struct HaplotypeLimits { unsigned target = 128, holdout = 2048, overflow = 8192; } haplotype_limits;
        unsigned max_holdout_depth = 2;
        Haplotype::MappingDomain::Size min_flank_pad = 30;
        boost::optional<Haplotype::NucleotideSequence::size_type> max_indicator_join_distance = boost::none;
        boost::optional<GenomicRegion::Distance> max_allele_distance = boost::none;
    };
    
    enum class Mode { allele, haplotype, allele_and_haplotype };
    
    class HaplotypeOverflow;
    class Builder;
    
    using HaplotypeBlock = HaplotypeTree::HaplotypeBlock;
    struct HaplotypePacket
    {
        HaplotypeBlock haplotypes;
        boost::optional<GenomicRegion> active_region, backtrack_region;
    };
    
    HaplotypeGenerator() = delete;
    
    HaplotypeGenerator(const ReferenceGenome& reference,
                       const MappableFlatSet<Variant>& candidates,
                       const ReadMap& reads,
                       boost::optional<const TemplateMap&> read_templates,
                       Policies policies);
    
    HaplotypeGenerator(const HaplotypeGenerator&)            = default;
    HaplotypeGenerator& operator=(const HaplotypeGenerator&) = default;
    HaplotypeGenerator(HaplotypeGenerator&&)                 = default;
    HaplotypeGenerator& operator=(HaplotypeGenerator&&)      = default;
    
    ~HaplotypeGenerator() = default;
    
    // Generates the next HaplotypePacket, this will always move the generated
    // active region forward. The next generated region is automatically calculated
    // (i.e the region given by peek_next_active_region) unless jump has
    // been called.
    HaplotypePacket generate();
    
    // Returns the next active region assuming jump is not called,
    // unless the generator is in holdout mode in which case boost::none is returned.
    boost::optional<GenomicRegion> peek_next_active_region() const;
    
    // Clears the history from the current generator, this has the effect of ensuring
    // the next haplotype set will not be lagged.
    void clear_progress() noexcept;
    
    // Unconditionally moves the next active region to within a specific interval.
    // Note this does not mean the next active region will be equal to the given
    // arguement, only that it will be contained by it.
    void jump(GenomicRegion to);
    
    // Discards the given set of haplotypes from generation, so that the 
    // next returned haplotype set will not contain these haplotypes as
    // sub-haplotypes.
    template <typename Container> void remove(const Container& haplotypes);
    
    // Remove all active haplotypes in the given region
    void remove(const GenomicRegion& region);
    
    // Returns true if calling remove will change the next active region.
    bool removal_has_impact() const;
    
    // Returns the number of haplotypes that can be removed from the generator
    // that will have a removal impact.
    unsigned max_removal_impact() const;
    
    // Discards any equivilant haplotypes that are not in the given set of
    // haplotypes.
    void collapse(const HaplotypeBlock& haplotypes);
    
    bool done() const noexcept;
    
    void add_lagging_exclusion_zone(GenomicRegion region);
    void clear_lagging_exclusion_zones() noexcept;
    
private:
    struct HoldoutSet
    {
        std::vector<Allele> alleles;
        GenomicRegion region;
        
        template <typename InputIt>
        HoldoutSet(InputIt first, InputIt last, GenomicRegion region)
        : alleles {first, last}
        , region {std::move(region)}
        {}
    };
    
    Policies policies_;
    
    HaplotypeTree tree_;
    GenomeWalker default_walker_, holdout_walker_;
    boost::optional<GenomeWalker> lagged_walker_;
    
    MappableFlatSet<Allele> alleles_;
    std::reference_wrapper<const ReadMap> reads_;
    boost::optional<const TemplateMap&> read_templates_;
    
    MappableFlatSet<GenomicRegion> lagging_exclusion_zones_;
    
    GenomicRegion active_region_;
    mutable boost::optional<GenomicRegion> next_active_region_;
    
    std::stack<HoldoutSet> active_holdouts_;
    boost::optional<GenomicRegion> holdout_region_;
    std::set<std::vector<ContigRegion>> previous_holdout_regions_;
    
    Allele rightmost_allele_;
    
    std::deque<HaplotypeBlock> haplotype_blocks_;
    std::vector<HaplotypeBlock> active_haplotype_blocks_;
    
    mutable boost::optional<logging::DebugLogger> debug_log_;
    mutable boost::optional<logging::TraceLogger> trace_log_;
    
    bool is_lagging_enabled() const noexcept;
    bool is_lagging_enabled(const GenomicRegion& region) const;
    bool is_backtracking_enabled() const noexcept;
    bool is_backtrack_active() const noexcept;
    bool is_active_region_lagged() const;
    void reset_next_active_region() const noexcept;
    GenomicRegion walk_from_active_region(const GenomeWalker& walker) const;
    GenomicRegion find_max_lagged_region() const;
    void update_next_active_region() const;
    void update_lagged_next_active_region() const;
    void prepare_for_next_active_region();
    void safe_clear_tree(const GenomicRegion& region);
    void populate_tree();
    void populate_tree_with_novel_alleles();
    void populate_tree_with_holdouts();
    unsigned extend_tree_with_cached_haplotypes(unsigned max_haplotypes);
    void extend_tree_with_cached_haplotypes();
    void clear_and_populate_tree_with_cached_haplotypes();
    bool can_populate_tree_with_cached_haplotype(unsigned max_haplotypes) const;
    bool can_populate_tree_with_cached_haplotype() const;
    bool should_populate_tree_with_cached_haplotype() const;
    std::size_t max_cached_haplotypes() const;
    bool in_holdout_mode() const noexcept;
    const GenomicRegion& top_holdout_region() const;
    bool can_try_extracting_holdouts(const GenomicRegion& region) const noexcept;
    bool try_extract_holdouts(GenomicRegion region);
    boost::optional<HoldoutSet> propose_new_holdout_set(GenomicRegion region) const;
    bool can_reintroduce_holdouts() const noexcept;
    void reintroduce_holdouts();
    void clear_holdouts() noexcept;
    void resolve_sandwich_insertion();
    GenomicRegion calculate_haplotype_region() const;
    bool is_finished() const noexcept;
    void cleanup_tree();
    void cache_active_haplotypes(const GenomicRegion& region);
    void cache_active_haplotypes();
    boost::optional<GenomicRegion> haplotype_block_region() const;
    boost::optional<GenomicRegion> backtrack_region() const;
};

class HaplotypeGenerator::HaplotypeOverflow : public std::runtime_error
{
public:
    HaplotypeOverflow(GenomicRegion region, std::size_t size);
    
    virtual ~HaplotypeOverflow() noexcept = default;
    
    virtual const char* what() const noexcept override;
    const GenomicRegion& region() const noexcept;
    std::size_t size() const noexcept;
    
private:
    GenomicRegion region_;
    std::size_t size_;
    std::string message_;
};

template <typename Container>
void HaplotypeGenerator::remove(const Container& haplotypes)
{
    if (done() || haplotypes.empty()) return;
    reset_next_active_region();
    prune_all(haplotypes, tree_);
    if (debug_log_)
        stream(*debug_log_) << "There are " << tree_.num_haplotypes() << " left in the tree after pruning "
                            << haplotypes.size() << " haplotypes";
}

class HaplotypeGenerator::Builder
{
public:
    using Policies = HaplotypeGenerator::Policies;
    
    Builder() = default;
    
    Builder(const Builder&)            = default;
    Builder& operator=(const Builder&) = default;
    Builder(Builder&&)                 = default;
    Builder& operator=(Builder&&)      = default;
    
    ~Builder() = default;
    
    Builder& set_lagging_policy(Policies::Lagging policy) noexcept;
    Builder& set_extension_policy(Policies::Extension policy) noexcept;
    Builder& set_backtrack_policy(Policies::Backtrack policy) noexcept;
    Builder& set_target_limit(unsigned n) noexcept;
    Builder& set_holdout_limit(unsigned n) noexcept;
    Builder& set_overflow_limit(unsigned n) noexcept;
    Builder& set_max_holdout_depth(unsigned n) noexcept;
    Builder& set_min_flank_pad(Haplotype::MappingDomain::Size n) noexcept;
    Builder& set_max_indicator_join_distance(Haplotype::NucleotideSequence::size_type n) noexcept;
    Builder& set_max_allele_distance(GenomicRegion::Distance gap) noexcept;
    
    HaplotypeGenerator
    build(const ReferenceGenome& reference,
          const MappableFlatSet<Variant>& candidates,
          const ReadMap& reads,
          boost::optional<const TemplateMap&> read_templates = boost::none) const;
    
private:
    Policies policies_;
};

} // namespace coretools

using coretools::HaplotypeGenerator;

} // namespace octopus

#endif
