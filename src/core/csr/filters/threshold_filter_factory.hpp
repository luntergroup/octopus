// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef threshold_filter_factory_hpp
#define threshold_filter_factory_hpp

#include <memory>
#include <vector>
#include <string>

#include <boost/optional.hpp>

#include "logging/progress_meter.hpp"
#include "variant_call_filter_factory.hpp"
#include "variant_call_filter.hpp"
#include "threshold_filter.hpp"
#include "somatic_threshold_filter.hpp"

namespace octopus { namespace csr {

class FacetFactory;

class ThresholdFilterFactory : public VariantCallFilterFactory
{
public:
    ThresholdFilterFactory() = default;
    
    ThresholdFilterFactory(std::string soft_expression);
    ThresholdFilterFactory(std::string hard_expression, std::string soft_expression);
    ThresholdFilterFactory(std::string germline_hard_expression, std::string germline_soft_expression,
                           std::string somatic_hard_expression, std::string somatic_soft_expression,
                           std::string refcall_hard_expression, std::string refcall_soft_expression);
    
    ThresholdFilterFactory(const ThresholdFilterFactory&)            = default;
    ThresholdFilterFactory& operator=(const ThresholdFilterFactory&) = default;
    ThresholdFilterFactory(ThresholdFilterFactory&&)                 = default;
    ThresholdFilterFactory& operator=(ThresholdFilterFactory&&)      = default;
    
    ~ThresholdFilterFactory() = default;

private:
    using Condition = ThresholdVariantCallFilter::Condition;
    using ConditionVectorPair = ThresholdVariantCallFilter::ConditionVectorPair;
    
    ConditionVectorPair germline_, somatic_, reference_;
    
    std::unique_ptr<VariantCallFilterFactory> do_clone() const override;
    std::unique_ptr<VariantCallFilter> do_make(FacetFactory facet_factory,
                                               VariantCallFilter::OutputOptions output_config,
                                               boost::optional<ProgressMeter&> progress,
                                               VariantCallFilter::ConcurrencyPolicy threading) const override;
};

} // namespace csr

using csr::ThresholdFilterFactory;

} // namespace octopus

#endif
