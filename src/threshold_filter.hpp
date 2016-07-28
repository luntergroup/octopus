//
//  threshold_filter.hpp
//  Octopus
//
//  Created by Daniel Cooke on 19/07/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef threshold_filter_hpp
#define threshold_filter_hpp

#include "variant_call_filter.hpp"

class VcfHeader;

namespace octopus { namespace CallFiltering
{
    class ThresholdVariantCallFilter : public VariantCallFilter
    {
    public:
        ThresholdVariantCallFilter() = delete;
        
        ThresholdVariantCallFilter(const ReferenceGenome& reference,
                                   const ReadPipe& read_pipe,
                                   std::vector<MeasureWrapper> measures,
                                   std::size_t max_read_buffer_size);
        
        ThresholdVariantCallFilter(const ThresholdVariantCallFilter&)            = delete;
        ThresholdVariantCallFilter& operator=(const ThresholdVariantCallFilter&) = delete;
        ThresholdVariantCallFilter(ThresholdVariantCallFilter&&)                 = default;
        ThresholdVariantCallFilter& operator=(ThresholdVariantCallFilter&&)      = default;
        
        virtual ~ThresholdVariantCallFilter() = default;
        
    private:
        virtual void annotate(VcfHeader& dest) const override;
        virtual Classification classify(const MeasureVector& call_measures) const override;
    };
} // namespace CallFiltering
} // namespace octopus

#endif /* threshold_filter_hpp */
