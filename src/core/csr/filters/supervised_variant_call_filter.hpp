//
//  supervised_variant_call_filter.hpp
//  octopus
//
//  Created by Daniel Cooke on 23/07/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef supervised_variant_call_filter_hpp
#define supervised_variant_call_filter_hpp

#include <deque>

#include "variant_call_filter.hpp"

namespace octopus { namespace csr
{
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

#endif /* supervised_variant_call_filter_hpp */
