//
//  genome_walker.hpp
//  Octopus
//
//  Created by Daniel Cooke on 17/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef genome_walker_hpp
#define genome_walker_hpp

#include <functional>

#include "common.hpp"
#include "mappable.hpp"
#include "allele.hpp"
#include "mappable_flat_set.hpp"

class GenomicRegion;

namespace octopus {

class Variant;
    
namespace coretools {

class GenomeWalker
{
public:
    using ContigName = GenomicRegion::ContigName;
    using Candidates = MappableFlatSet<Allele>;
    
    enum class IndicatorPolicy
    {
        IncludeNone,
        IncludeIfLinkableToNovelRegion,
        IncludeIfSharedWithNovelRegion,
        IncludeAll
    };
    
    enum class ExtensionPolicy
    {
        IncludeIfWithinReadLengthOfFirstIncluded,
        IncludeIfSharedWithFrontier,
        NoLimit
    };
    
    GenomeWalker() = delete;
    
    GenomeWalker(unsigned max_included,
                 IndicatorPolicy indicator_policy = IndicatorPolicy::IncludeNone,
                 ExtensionPolicy extension_policy = ExtensionPolicy::IncludeIfSharedWithFrontier);
    
    GenomeWalker(const GenomeWalker&)            = default;
    GenomeWalker& operator=(const GenomeWalker&) = default;
    GenomeWalker(GenomeWalker&&)                 = default;
    GenomeWalker& operator=(GenomeWalker&&)      = default;
    
    ~GenomeWalker() = default;
    
    GenomicRegion walk(const ContigName& contig, const ReadMap& reads,
                       const Candidates& candidates) const;
    
    GenomicRegion walk(const GenomicRegion& previous_region, const ReadMap& reads,
                       const Candidates& candidates) const;
    
private:
    unsigned max_included_;
    
    IndicatorPolicy indicator_policy_;
    ExtensionPolicy extension_policy_;
    
    using CandidateIterator = Candidates::const_iterator;
};

} // namespace coretools
} // namespace octopus

#endif /* genome_walker_hpp */
