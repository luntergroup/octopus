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
#include "mappable_flat_multi_set.hpp"
#include "allele.hpp"

#include "mappable_flat_set.hpp"

class GenomicRegion;
class Variant;

namespace Octopus
{
class GenomeWalker
{
public:
    using ContigNameType = GenomicRegion::ContigNameType;
    using Candidates     = MappableFlatSet<Allele>;
    
    enum class IndicatorLimit { SharedWithPreviousRegion, NoLimit };
    enum class ExtensionLimit { WithinReadLengthOfFirstIncluded, SharedWithFrontier, NoLimit };
    
    GenomeWalker() = delete;
    
    explicit GenomeWalker(unsigned max_indicators,
                          unsigned max_included,
                          IndicatorLimit indicator_limit = IndicatorLimit::SharedWithPreviousRegion,
                          ExtensionLimit extension_limit = ExtensionLimit::NoLimit);
    
    ~GenomeWalker() = default;
    
    GenomeWalker(const GenomeWalker&)            = default;
    GenomeWalker& operator=(const GenomeWalker&) = default;
    GenomeWalker(GenomeWalker&&)                 = default;
    GenomeWalker& operator=(GenomeWalker&&)      = default;
    
    GenomicRegion
    walk(const ContigNameType& contig, const ReadMap& reads, const Candidates& candidates) const;
    
    GenomicRegion
    walk(const GenomicRegion& previous_region, const ReadMap& reads, const Candidates& candidates) const;
    
private:
    const unsigned max_indicators_;
    const unsigned max_included_;
    
    const IndicatorLimit indicator_limit_;
    const ExtensionLimit extension_limit_;
    
    using CandidateIterator = Candidates::const_iterator;
};
} // namespace Octopus

#endif /* genome_walker_hpp */
