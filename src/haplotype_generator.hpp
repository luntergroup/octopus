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
#include <unordered_map>

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
        HaplotypeGenerator() = delete;
        
        explicit HaplotypeGenerator(const GenomicRegion& window, const ReferenceGenome& reference,
                                    const MappableFlatSet<Variant>& candidates, const ReadMap& reads,
                                    unsigned max_haplotypes, bool allow_lagging);
        
        ~HaplotypeGenerator() = default;
        
        HaplotypeGenerator(const HaplotypeGenerator&)            = default;
        HaplotypeGenerator& operator=(const HaplotypeGenerator&) = default;
        HaplotypeGenerator(HaplotypeGenerator&&)                 = default;
        HaplotypeGenerator& operator=(HaplotypeGenerator&&)      = default;
        
        const GenomicRegion& tell_next_active_region() const;
        
        std::pair<std::vector<Haplotype>, GenomicRegion> progress();
        
        void clear_progress() noexcept;
        
        void uniquely_keep(const std::vector<Haplotype>& haplotypes);
        void remove(const std::vector<Haplotype>& haplotypes);
        void remove(const std::vector<std::reference_wrapper<const Haplotype>>& haplotypes);
        
        void force_forward(GenomicRegion to);
        
    private:
        HaplotypeTree tree_;
        GenomeWalker walker_;
        boost::optional<GenomeWalker> lagged_walker_;
        
        MappableFlatSet<Allele> alleles_;
        std::reference_wrapper<const ReadMap> reads_;
        
        GenomicRegion current_active_region_;
        mutable boost::optional<GenomicRegion> next_active_region_;
        
        unsigned max_haplotypes_, hard_max_haplotypes_ = 150'000;
        
        MappableFlatMultiSet<Allele> holdout_set_;
        
        boost::optional<Allele> rightmost_allele_;
        
        bool is_lagged() const noexcept;
        bool is_active_region_lagged() const;
        
        void update_next_active_region() const;
        MappableFlatMultiSet<Allele> compute_holdout_set(const GenomicRegion& active_region) const;
        GenomicRegion calculate_haplotype_region() const;
    };
} // namespace Octopus

#endif /* haplotype_generator_hpp */
