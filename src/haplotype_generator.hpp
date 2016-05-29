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
    
    class HaplotypeOverflowError;
    
    class Builder;
    
    HaplotypeGenerator() = delete;
    
    HaplotypeGenerator(const GenomicRegion& window, const ReferenceGenome& reference,
                       const MappableFlatSet<Variant>& candidates, const ReadMap& reads,
                       unsigned soft_max_haplotypes, unsigned hard_max_haplotypes,
                       LaggingPolicy lagging = LaggingPolicy::None);
    
    ~HaplotypeGenerator() = default;
    
    HaplotypeGenerator(const HaplotypeGenerator&)            = default;
    HaplotypeGenerator& operator=(const HaplotypeGenerator&) = default;
    HaplotypeGenerator(HaplotypeGenerator&&)                 = default;
    HaplotypeGenerator& operator=(HaplotypeGenerator&&)      = default;
    
    GenomicRegion tell_next_active_region() const;
    
    void progress(GenomicRegion to);
    
    void stop() noexcept;
    
    std::pair<std::vector<Haplotype>, GenomicRegion> generate();
    
    bool removal_has_impact() const;
    
    template <typename Container> void remove_duplicates(const Container& haplotypes);
    template <typename Container> void remove(const Container& haplotypes);
    
private:
    HaplotypeTree tree_;
    GenomeWalker walker_;
    boost::optional<GenomeWalker> lagged_walker_;
    
    unsigned soft_max_haplotypes_, hard_max_haplotypes_;
    
    MappableFlatSet<Allele> alleles_;
    std::reference_wrapper<const ReadMap> reads_;
    
    GenomicRegion current_active_region_;
    
    mutable boost::optional<GenomicRegion> next_active_region_;
    
    mutable MappableFlatMultiSet<Allele> holdout_set_;
    mutable boost::optional<GenomicRegion> current_holdout_region_;
    mutable std::unordered_set<GenomicRegion> previous_holdout_regions_;
    
    boost::optional<Allele> rightmost_allele_;
    
    bool is_lagged() const noexcept;
    bool is_active_region_lagged() const;
    
    void reset_next_active_region() const noexcept;
    void update_next_active_region() const;
    
    void set_holdout_set(const GenomicRegion& active_region);
    void rientroduce_holdout_set();
    
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
        alleles_.erase_overlapped(current_active_region_);
    } else {
        prune_all(haplotypes, tree_);
    }
}

class HaplotypeGenerator::Builder
{
public:
    using LaggingPolicy = HaplotypeGenerator::LaggingPolicy;
    
    Builder();
    
    ~Builder() = default;
    
    Builder(const Builder&)            = default;
    Builder& operator=(const Builder&) = default;
    Builder(Builder&&)                 = default;
    Builder& operator=(Builder&&)      = default;
    
    Builder& set_soft_max_haplotypes(unsigned n) noexcept;
    Builder& set_soft_hard_haplotypes(unsigned n) noexcept;
    Builder& set_lagging_policy(LaggingPolicy policy) noexcept;
    
    HaplotypeGenerator build(const ReferenceGenome& reference, const GenomicRegion& window,
                             const MappableFlatSet<Variant>& candidates,
                             const ReadMap& reads) const;
    
private:
    unsigned soft_max_haplotypes_, hard_max_haplotypes_;
    LaggingPolicy lagging_policy_;
};
} // namespace Octopus

#endif /* haplotype_generator_hpp */
