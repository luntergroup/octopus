//
//  variant_call_filter.hpp
//  Octopus
//
//  Created by Daniel Cooke on 31/05/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef variant_call_filter_hpp
#define variant_call_filter_hpp

#include <vector>

#include "vcf_record.hpp"

namespace Octopus
{
    class VariantCallFilter
    {
    public:
        VariantCallFilter() = delete;
        
        virtual ~VariantCallFilter() = default;
        
        VariantCallFilter(const VariantCallFilter&)            = delete;
        VariantCallFilter& operator=(const VariantCallFilter&) = delete;
        VariantCallFilter(VariantCallFilter&&)                 = default;
        VariantCallFilter& operator=(VariantCallFilter&&)      = default;
        
        virtual void apply_filters(std::vector<VcfRecord>& calls) const;
        
    private:
        
    };
} // namespace Octopus

#endif /* variant_call_filter_hpp */
