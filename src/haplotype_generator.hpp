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
#include <unordered_set>
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
    enum class LaggingPolicy { None, Conservative, Aggressive };
    
    struct HaplotypeLimits
    {
        unsigned soft_max, holdout_max, hard_max;
    };
    
    class HaplotypeOverflowError;
    
    class Builder;
    
    using HaplotypePacket = std::tuple<std::vector<Haplotype>, GenomicRegion, bool>;
    
    HaplotypeGenerator() = delete;
    
    HaplotypeGenerator(const GenomicRegion& window, const ReferenceGenome& reference,
                       const MappableFlatSet<Variant>& candidates, const ReadMap& reads,
                       HaplotypeLimits haplotype_limits,
                       LaggingPolicy lagging = LaggingPolicy::None,
                       Haplotype::SizeType min_pad = 30);
    
    ~HaplotypeGenerator() = default;
    
    HaplotypeGenerator(const HaplotypeGenerator&)            = default;
    HaplotypeGenerator& operator=(const HaplotypeGenerator&) = default;
    HaplotypeGenerator(HaplotypeGenerator&&)                 = default;
    HaplotypeGenerator& operator=(HaplotypeGenerator&&)      = default;
    
    GenomicRegion tell_next_active_region() const;
    
    void progress(GenomicRegion to);
    
    void stop() noexcept;
    
    HaplotypePacket generate();
    
    bool removal_has_impact() const;
    unsigned max_removal_impact() const;
    
    template <typename Container> void remove_duplicates(const Container& haplotypes);
    template <typename Container> void remove(const Container& haplotypes);
    
private:
    HaplotypeTree tree_;
    GenomeWalker walker_;
    boost::optional<GenomeWalker> lagged_walker_;
    
    HaplotypeLimits haplotype_limits_;
    
    Haplotype::SizeType min_pad_ = 30;
    
    MappableFlatSet<Allele> alleles_;
    std::reference_wrapper<const ReadMap> reads_;
    
    GenomicRegion current_active_region_;
    
    mutable boost::optional<GenomicRegion> next_active_region_;
    
    mutable MappableFlatSet<Allele> holdout_set_;
    mutable boost::optional<GenomicRegion> current_holdout_region_;
    mutable std::unordered_set<GenomicRegion> previous_holdout_regions_;
    
    boost::optional<Allele> rightmost_allele_;
    
    bool is_lagging_enabled() const noexcept;
    bool is_active_region_lagged() const;
    
    void reset_next_active_region() const noexcept;
    void update_next_active_region() const;
    
    void set_holdout_set(const GenomicRegion& active_region);
    bool try_reintroducing_holdout_set();
    
    GenomicRegion calculate_haplotype_region() const;
};

class HaplotypeGenerator::HaplotypeOverflowError : public std::runtime_error
{
public:
    HaplotypeOverflowError(GenomicRegion region);
    
    virtual const char* what() const noexcept override;
    
    const GenomicRegion& region() const noexcept;
    
    unsigned overflow_size() const noexcept;
    
private:
    GenomicRegion region_;
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
        if (holdout_set_.empty()) {
            alleles_.erase_overlapped(current_active_region_);
        } else {
            // TODO: in this case we must be more selective and only erase those alleles
            // which are not present in the remaining haplotype set
            alleles_.erase_overlapped(current_active_region_);
        }
    } else {
        prune_all(haplotypes, tree_);
    }
}

class HaplotypeGenerator::Builder
{
public:
    using LaggingPolicy = HaplotypeGenerator::LaggingPolicy;
    
    Builder() = default;
    
    ~Builder() = default;
    
    Builder(const Builder&)            = default;
    Builder& operator=(const Builder&) = default;
    Builder(Builder&&)                 = default;
    Builder& operator=(Builder&&)      = default;
    
    Builder& set_soft_max_haplotypes(unsigned n) noexcept;
    Builder& set_hard_max_haplotypes(unsigned n) noexcept;
    Builder& set_lagging_policy(LaggingPolicy policy) noexcept;
    Builder& set_min_pad(Haplotype::SizeType val) noexcept;
    
    HaplotypeGenerator build(const ReferenceGenome& reference, const GenomicRegion& window,
                             const MappableFlatSet<Variant>& candidates, const ReadMap& reads) const;
    
private:
    HaplotypeGenerator::HaplotypeLimits haplotype_limits_ = {128, 2048, 16384};
    LaggingPolicy lagging_policy_ = LaggingPolicy::None;
    Haplotype::SizeType min_pad_ = 30;
};
} // namespace Octopus

#endif /* haplotype_generator_hpp */
