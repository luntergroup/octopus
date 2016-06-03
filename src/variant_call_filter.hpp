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

#include "read_manager.hpp"
#include "vcf_reader.hpp"
#include "vcf_writer.hpp"

namespace Octopus
{
    class VariantCallFilter
    {
    public:
        VariantCallFilter() = delete;
        
        VariantCallFilter(ReadManager& read_manager);
        
        virtual ~VariantCallFilter() = default;
        
        VariantCallFilter(const VariantCallFilter&)            = delete;
        VariantCallFilter& operator=(const VariantCallFilter&) = delete;
        VariantCallFilter(VariantCallFilter&&)                 = default;
        VariantCallFilter& operator=(VariantCallFilter&&)      = default;
        
        void filter(VcfReader& source, VcfWriter& dest);
        
    private:
        std::reference_wrapper<ReadManager> read_manager_;
    };
} // namespace Octopus

#endif /* variant_call_filter_hpp */
