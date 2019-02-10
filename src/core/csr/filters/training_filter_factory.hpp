// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef training_filter_factory_hpp
#define training_filter_factory_hpp

#include <memory>
#include <vector>
#include <string>
#include <set>

#include <boost/optional.hpp>

#include "logging/progress_meter.hpp"
#include "variant_call_filter_factory.hpp"
#include "variant_call_filter.hpp"
#include "../measures/measure.hpp"

namespace octopus { namespace csr {

class FacetFactory;

class TrainingFilterFactory : public VariantCallFilterFactory
{
public:
    TrainingFilterFactory() = default;
    
    TrainingFilterFactory(const std::set<std::string>& measures);
    
    TrainingFilterFactory(const TrainingFilterFactory&)            = default;
    TrainingFilterFactory& operator=(const TrainingFilterFactory&) = default;
    TrainingFilterFactory(TrainingFilterFactory&&)                 = default;
    TrainingFilterFactory& operator=(TrainingFilterFactory&&)      = default;
    
    ~TrainingFilterFactory() = default;

private:
    std::vector<MeasureWrapper> measures_;
    
    std::unique_ptr<VariantCallFilterFactory> do_clone() const override;
    std::unique_ptr<VariantCallFilter> do_make(FacetFactory facet_factory,
                                               VariantCallFilter::OutputOptions output_config,
                                               boost::optional<ProgressMeter&> progress,
                                               VariantCallFilter::ConcurrencyPolicy threading) const override;
};

} // namespace csr

using csr::TrainingFilterFactory;

} // namespace octopus

#endif
