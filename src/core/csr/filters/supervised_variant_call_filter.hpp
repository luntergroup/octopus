// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef supervised_variant_call_filter_hpp
#define supervised_variant_call_filter_hpp

#include <deque>

#include "variant_call_filter.hpp"
#include "io/reference/reference_genome.hpp"
#include "readpipe/read_pipe_fwd.hpp"

namespace octopus { namespace csr  {

class SupervisedVariantCallFilter : public VariantCallFilter
{
public:
    SupervisedVariantCallFilter() = delete;
    
    SupervisedVariantCallFilter(const ReferenceGenome& reference,
                                const ReadPipe& read_pipe,
                                std::vector<MeasureWrapper> measures,
                                std::size_t max_read_buffer_size);
    
    SupervisedVariantCallFilter(const SupervisedVariantCallFilter&)            = delete;
    SupervisedVariantCallFilter& operator=(const SupervisedVariantCallFilter&) = delete;
    SupervisedVariantCallFilter(SupervisedVariantCallFilter&&)                 = default;
    SupervisedVariantCallFilter& operator=(SupervisedVariantCallFilter&&)      = default;
    
    virtual ~SupervisedVariantCallFilter() = default;
    
protected:
    struct TrainingCall
    {
        Variant variant;
        MeasureVector measures;
        double confidence;
    };
    
    std::deque<TrainingCall> training_calls_;
    
private:
    virtual bool is_online() const noexcept = 0;
};

} // namespace csr
} // namespace octopus

#endif
