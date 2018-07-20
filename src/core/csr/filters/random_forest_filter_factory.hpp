// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef random_forest_filter_factory_hpp
#define random_forest_filter_factory_hpp

#include <memory>
#include <vector>
#include <string>

#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

#include "basics/phred.hpp"
#include "logging/progress_meter.hpp"
#include "../measures/measure.hpp"
#include "variant_call_filter_factory.hpp"
#include "variant_call_filter.hpp"
#include "random_forest_filter.hpp"

namespace octopus { namespace csr {

class FacetFactory;

class RandomForestFilterFactory : public VariantCallFilterFactory
{
public:
    using Path = RandomForestFilter::Path;
    enum class ForestType { germline, somatic, denovo };
    
    RandomForestFilterFactory();
    RandomForestFilterFactory(Path ranger_forest, Path temp_directory, ForestType type = ForestType::germline);
    RandomForestFilterFactory(Path germline_ranger_forest, Path somatic_ranger_forest, Path temp_directory);
    
    RandomForestFilterFactory(const RandomForestFilterFactory&)            = default;
    RandomForestFilterFactory& operator=(const RandomForestFilterFactory&) = default;
    RandomForestFilterFactory(RandomForestFilterFactory&&)                 = default;
    RandomForestFilterFactory& operator=(RandomForestFilterFactory&&)      = default;
    
    ~RandomForestFilterFactory() = default;
    
    std::vector<MeasureWrapper> measures() const;
    
    void set_min_forest_quality(Phred<double> quality);

private:
    std::vector<MeasureWrapper> measures_;
    std::vector<Path> ranger_forests_;
    std::vector<ForestType> forest_types_;
    Path temp_directory_;
    boost::optional<Phred<double>> min_forest_quality_;
    
    std::unique_ptr<VariantCallFilterFactory> do_clone() const override;
    std::unique_ptr<VariantCallFilter> do_make(FacetFactory facet_factory,
                                               VariantCallFilter::OutputOptions output_config,
                                               boost::optional<ProgressMeter&> progress,
                                               VariantCallFilter::ConcurrencyPolicy threading) const override;
};

} // namespace csr

using csr::RandomForestFilterFactory;

} // namespace octopus

#endif
