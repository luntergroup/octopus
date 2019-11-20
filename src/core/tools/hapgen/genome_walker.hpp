// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef genome_walker_hpp
#define genome_walker_hpp

#include <functional>

#include "config/common.hpp"
#include "concepts/mappable.hpp"
#include "core/types/allele.hpp"
#include "containers/mappable_flat_set.hpp"

namespace octopus {

class GenomicRegion;
class Variant;

namespace coretools {

class GenomeWalker
{
public:
    using AlleleSet = MappableFlatSet<Allele>;
    
    enum class IndicatorPolicy
    {
        includeNone,
        includeIfSharedWithNovelRegion,
        includeIfLinkableToNovelRegion,
        includeAll
    };
    enum class ExtensionPolicy
    {
        includeIfWithinReadLengthOfFirstIncluded,
        includeIfAllSamplesSharedWithFrontier,
        includeIfAnySampleSharedWithFrontier,
        noLimit
    };
    
    enum class ReadTemplatePolicy
    {
        none,
        indicators,
        extension,
        indicators_and_extension
    };
    
    GenomeWalker() = delete;
    
    GenomeWalker(unsigned max_included,
                 IndicatorPolicy indicator_policy = IndicatorPolicy::includeNone,
                 ExtensionPolicy extension_policy = ExtensionPolicy::includeIfAnySampleSharedWithFrontier,
                 ReadTemplatePolicy read_template_policy = ReadTemplatePolicy::indicators_and_extension);
    
    GenomeWalker(const GenomeWalker&)            = default;
    GenomeWalker& operator=(const GenomeWalker&) = default;
    GenomeWalker(GenomeWalker&&)                 = default;
    GenomeWalker& operator=(GenomeWalker&&)      = default;
    
    ~GenomeWalker() = default;
    
    GenomicRegion
    walk(const GenomicRegion& previous_region,
         const ReadMap& reads,
         const AlleleSet& alleles,
         boost::optional<const TemplateMap&> read_templates = boost::none) const;
    
private:
    unsigned max_included_;
    IndicatorPolicy indicator_policy_;
    ExtensionPolicy extension_policy_;
    ReadTemplatePolicy read_template_policy_;
    
    bool can_extend(const Allele& active, const Allele& novel,
                    const ReadMap& reads, boost::optional<const TemplateMap&> read_templates) const;
};

} // namespace coretools
} // namespace octopus

#endif
