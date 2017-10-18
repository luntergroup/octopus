// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef threshold_filter_hpp
#define threshold_filter_hpp

#include "variant_call_filter.hpp"
#include "io/reference/reference_genome.hpp"
#include "readpipe/read_pipe_fwd.hpp"

namespace octopus {

class VcfHeader;

namespace csr {

class ThresholdVariantCallFilter : public VariantCallFilter
{
public:
    struct Threshold
    {
        virtual bool operator()(Measure::ResultType value) const noexcept { return true; }
        virtual ~Threshold() = default;
    };
    
    ThresholdVariantCallFilter() = delete;
    
    ThresholdVariantCallFilter(const ReferenceGenome& reference,
                               const ReadPipe& read_pipe,
                               std::vector<MeasureWrapper> measures,
                               std::vector<std::unique_ptr<Threshold>> thresholds);
    
    ThresholdVariantCallFilter(const ThresholdVariantCallFilter&)            = delete;
    ThresholdVariantCallFilter& operator=(const ThresholdVariantCallFilter&) = delete;
    ThresholdVariantCallFilter(ThresholdVariantCallFilter&&)                 = default;
    ThresholdVariantCallFilter& operator=(ThresholdVariantCallFilter&&)      = default;
    
    virtual ~ThresholdVariantCallFilter() = default;

private:
    std::reference_wrapper<const ReadPipe> read_pipe_;
    std::vector<std::unique_ptr<Threshold>> thresholds_;
    
    virtual void annotate(VcfHeader& dest) const override;
    virtual Classification classify(const MeasureVector& measures) const override;
    
    bool passes_all_filters(const MeasureVector& measures) const;
};

} // namespace csr
} // namespace octopus

#endif /* threshold_filter_hpp */
