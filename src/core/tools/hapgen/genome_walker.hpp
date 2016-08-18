// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef genome_walker_hpp
#define genome_walker_hpp

#include <functional>

#include <config/common.hpp>
#include <concepts/mappable.hpp>
#include <core/types/allele.hpp>
#include <containers/mappable_flat_set.hpp>

namespace octopus {

class GenomicRegion;
class Variant;
    
namespace coretools {

class GenomeWalker
{
public:
    using ContigName = GenomicRegion::ContigName;
    using Candidates = MappableFlatSet<Allele>;
    
    enum class IndicatorPolicy
    {
        includeNone,
        includeIfLinkableToNovelRegion,
        includeIfSharedWithNovelRegion,
        includeAll
    };
    
    enum class ExtensionPolicy
    {
        includeIfWithinReadLengthOfFirstIncluded,
        includeIfSharedWithFrontier,
        noLimit
    };
    
    GenomeWalker() = delete;
    
    GenomeWalker(unsigned max_included,
                 IndicatorPolicy indicator_policy = IndicatorPolicy::includeNone,
                 ExtensionPolicy extension_policy = ExtensionPolicy::includeIfSharedWithFrontier);
    
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

#endif
