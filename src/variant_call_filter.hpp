//
//  variant_call_filter.hpp
//  Octopus
//
//  Created by Daniel Cooke on 31/05/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef variant_call_filter_hpp
#define variant_call_filter_hpp

#include <functional>
#include <map>

#include "common.hpp"
#include "reference_genome.hpp"
#include "read_manager.hpp"

class GenomicRegion;
class VcfReader;
class VcfWriter;

namespace Octopus
{
    class VariantCallFilter
    {
    public:
        using ContigOrder = std::function<bool(const GenomicRegion&, const GenomicRegion&)>;
        using RegionMap   = std::map<ContigNameType, std::vector<GenomicRegion>, ContigOrder>;
        
        VariantCallFilter() = delete;
        
        VariantCallFilter(const ReferenceGenome& reference, const ReadManager& read_manager);
        
        virtual ~VariantCallFilter() = default;
        
        VariantCallFilter(const VariantCallFilter&)            = delete;
        VariantCallFilter& operator=(const VariantCallFilter&) = delete;
        VariantCallFilter(VariantCallFilter&&)                 = default;
        VariantCallFilter& operator=(VariantCallFilter&&)      = default;
        
        void filter(const VcfReader& source, VcfWriter& dest, const RegionMap& regions);
        
    private:
        std::reference_wrapper<const ReferenceGenome> reference_;
        std::reference_wrapper<const ReadManager> read_manager_;
    };
} // namespace Octopus

#endif /* variant_call_filter_hpp */
