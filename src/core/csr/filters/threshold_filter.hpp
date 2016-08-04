// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef threshold_filter_hpp
#define threshold_filter_hpp

#include "variant_call_filter.hpp"

namespace octopus {

class VcfHeader;

namespace csr
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
} // namespace csr
} // namespace octopus

#endif /* threshold_filter_hpp */
