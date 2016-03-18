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
                                    unsigned max_haplotypes, unsigned max_indicators);
        
        ~HaplotypeGenerator() = default;
        
        HaplotypeGenerator(const HaplotypeGenerator&)            = default;
        HaplotypeGenerator& operator=(const HaplotypeGenerator&) = default;
        HaplotypeGenerator(HaplotypeGenerator&&)                 = default;
        HaplotypeGenerator& operator=(HaplotypeGenerator&&)      = default;
        
        GenomicRegion tell_next_active_region() const;
        
        std::pair<std::vector<Haplotype>, GenomicRegion> progress();
        
        void keep_haplotypes(const std::vector<Haplotype>& haplotypes);
        void remove_haplotypes(const std::vector<Haplotype>& haplotypes);
        
        void force_forward(GenomicRegion to);
        
    private:
        HaplotypeTree tree_;
        GenomeWalker walker_;
        
        MappableFlatSet<Allele> alleles_;
        std::reference_wrapper<const ReadMap> reads_;
        
        GenomicRegion current_active_region_;
        mutable boost::optional<GenomicRegion> next_active_region_;
        
        unsigned max_haplotypes_, hard_max_haplotypes_ = 150'000;
        
        bool is_lagged_;
        
        MappableFlatMultiSet<Allele> holdout_set_;
        
        mutable std::unordered_map<Allele, unsigned> active_allele_counts_;
        
        MappableFlatMultiSet<Allele> compute_holdout_set(const GenomicRegion& active_region) const;
        GenomicRegion calculate_haplotype_region() const;
    };
} // namespace Octopus

#endif /* haplotype_generator_hpp */
