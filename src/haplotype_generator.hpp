//
//  haplotype_generator.hpp
//  Octopus
//
//  Created by Daniel Cooke on 20/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef haplotype_generator_hpp
#define haplotype_generator_hpp

#include <vector>
#include <tuple>
#include <functional>
#include <utility>
#include <stack>
#include <stdexcept>

#include <boost/optional.hpp>

#include "common.hpp"
#include "genome_walker.hpp"
#include "haplotype_tree.hpp"
#include "mappable_flat_set.hpp"
#include "mappable_flat_multi_set.hpp"
#include "genomic_region.hpp"
#include "allele.hpp"

class Variant;
class Haplotype;

namespace Octopus
{
class HaplotypeGenerator
{
public:
    struct Policies
    {
        enum class Lagging { None, Conservative, Aggressive } lagging = Lagging::None;
        
        struct HaplotypeLimits { unsigned target = 128, holdout = 2048, overflow = 8192; } haplotype_limits;
        
        unsigned max_holdout_depth = 2;
    };
    
    class HaplotypeOverflow;
    class Builder;
    
    using HaplotypePacket = std::pair<std::vector<Haplotype>, GenomicRegion>;
    
    HaplotypeGenerator() = delete;
    
    HaplotypeGenerator(const ReferenceGenome& reference,
                       const MappableFlatSet<Variant>& candidates,
                       const ReadMap& reads,
                       Policies policies,
                       Haplotype::SizeType min_flank_pad = 30);
    
    HaplotypeGenerator(const HaplotypeGenerator&)            = default;
    HaplotypeGenerator& operator=(const HaplotypeGenerator&) = default;
    HaplotypeGenerator(HaplotypeGenerator&&)                 = default;
    HaplotypeGenerator& operator=(HaplotypeGenerator&&)      = default;
    
    ~HaplotypeGenerator() = default;
    
    HaplotypePacket generate();
    
    boost::optional<GenomicRegion> peek_next_active_region() const;
    
    void force_progress(GenomicRegion to);
    
    void stop() noexcept;
    
    bool removal_has_impact() const;
    unsigned max_removal_impact() const;
    
    template <typename Container> void remove_duplicates(const Container& haplotypes);
    template <typename Container> void remove(const Container& haplotypes);
    
private:
    Policies policies_;
    Haplotype::SizeType min_flank_pad_;
    
    HaplotypeTree tree_;
    GenomeWalker default_walker_, holdout_walker_;
    boost::optional<GenomeWalker> lagged_walker_;
    
    MappableFlatSet<Allele> alleles_;
    std::reference_wrapper<const ReadMap> reads_;
    
    GenomicRegion active_region_;
    mutable boost::optional<GenomicRegion> next_active_region_;
    
    struct HoldoutSet
    {
        template <typename InputIt>
        HoldoutSet(InputIt first, InputIt last, GenomicRegion region)
        : alleles {first, last}, region {std::move(region)} {}
        std::vector<Allele> alleles;
        GenomicRegion region;
    };
    
    mutable std::stack<HoldoutSet> active_holdouts_;
    mutable boost::optional<GenomicRegion> holdout_region_;
    
    Allele rightmost_allele_;
    
    bool is_lagging_enabled() const noexcept;
    bool is_active_region_lagged() const;
    
    void reset_next_active_region() const noexcept;
    void update_next_active_region() const;
    
    void progress(GenomicRegion to);
    
    bool in_holdout_mode() const noexcept;
    bool can_extract_holdouts(const GenomicRegion& region) const noexcept;
    void extract_holdouts(GenomicRegion region);
    bool can_reintroduce_holdouts() const noexcept;
    void reintroduce_holdouts();
    void clear_holdouts() noexcept;
    
    GenomicRegion calculate_haplotype_region() const;
};

class HaplotypeGenerator::HaplotypeOverflow : public std::runtime_error
{
public:
    HaplotypeOverflow(GenomicRegion region, unsigned size);
    
    virtual ~HaplotypeOverflow() noexcept = default;
    
    virtual const char* what() const noexcept override;
    
    const GenomicRegion& region() const noexcept;
    
    unsigned size() const noexcept;
    
private:
    GenomicRegion region_;
    unsigned size_;
    std::string message_;
};

template <typename Container>
void HaplotypeGenerator::remove_duplicates(const Container& haplotypes)
{
    reset_next_active_region();
    
    if (!is_active_region_lagged() || tree_.num_haplotypes() == haplotypes.size()) {
        return;
    }
    
    prune_unique(haplotypes, tree_);
}

template <typename Container>
void HaplotypeGenerator::remove(const Container& haplotypes)
{
    if (haplotypes.empty()) return;
    
    reset_next_active_region();
    
    if (!is_active_region_lagged() || haplotypes.size() == tree_.num_haplotypes()) {
        tree_.clear();
        if (!in_holdout_mode()) {
            alleles_.erase_overlapped(active_region_);
        } else {
            // TODO: in this case we must be more selective and only erase those alleles
            // which are not present in the remaining haplotype set
            alleles_.erase_overlapped(active_region_);
        }
    } else {
        prune_all(haplotypes, tree_);
    }
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
    
    Builder& set_target_limit(unsigned n) noexcept;
    Builder& set_holdout_limit(unsigned n) noexcept;
    Builder& set_overflow_limit(unsigned n) noexcept;
    
    Builder& set_max_holdout_depth(unsigned n) noexcept;
    
    Builder& set_min_flank_pad(Haplotype::SizeType n) noexcept;
    
    HaplotypeGenerator build(const ReferenceGenome& reference,
                             const MappableFlatSet<Variant>& candidates,
                             const ReadMap& reads) const;
    
private:
    Policies policies_;
    Haplotype::SizeType min_flank_pad_ = 30;
};
} // namespace Octopus

#endif /* haplotype_generator_hpp */
