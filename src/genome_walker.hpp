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
#include "mappable_set.hpp"

class GenomicRegion;
class Variant;

namespace Octopus
{
class GenomeWalker
{
public:
    using ContigNameType = GenomicRegion::ContigNameType;
    using Candidates     = MappableSet<Variant>;
    
    enum class IndicatorLimit { SharedWithPreviousRegion, NoLimit };
    enum class ExtensionLimit { WithinReadLengthOfFirstIncluded, SharedWithFrontier, NoLimit };
    enum class ExpansionLimit { UpToExcluded, WithinReadLength, UpToExcludedWithinReadLength,
                                NoExpansion };
    
    GenomeWalker() = delete;
    
    explicit GenomeWalker(unsigned max_indicators,
                          unsigned max_included,
                          IndicatorLimit indicator_limit = IndicatorLimit::SharedWithPreviousRegion,
                          ExtensionLimit extension_limit = ExtensionLimit::SharedWithFrontier,
                          ExpansionLimit expansion_limit = ExpansionLimit::UpToExcludedWithinReadLength);
    
    ~GenomeWalker() = default;
    
    GenomeWalker(const GenomeWalker&)            = default;
    GenomeWalker& operator=(const GenomeWalker&) = default;
    GenomeWalker(GenomeWalker&&)                 = default;
    GenomeWalker& operator=(GenomeWalker&&)      = default;
    
    GenomicRegion walk(const ContigNameType& contig, const ReadMap& reads,
                       const Candidates& candidates);
    
    GenomicRegion walk(const GenomicRegion& previous_region, const ReadMap& reads,
                       const Candidates& candidates);
    
private:
    const unsigned max_indicators_;
    const unsigned max_included_;
    
    const IndicatorLimit indicator_limit_;
    const ExtensionLimit extension_limit_;
    const ExpansionLimit expansion_limit_;
    
    using CandidateIterator = Candidates::const_iterator;
    
    struct CandidateRanges
    {
        CandidateRanges() = delete;
        CandidateRanges(CandidateIterator first_previous_itr, CandidateIterator first_included_itr,
                        CandidateIterator first_excluded_itr, CandidateIterator last_itr);
        CandidateIterator first_previous_itr, first_included_itr, first_excluded_itr, last_itr;
    };
    
    std::function<GenomicRegion(CandidateRanges, const ReadMap&)> expander_;
};
} // namespace Octopus

#endif /* genome_walker_hpp */
