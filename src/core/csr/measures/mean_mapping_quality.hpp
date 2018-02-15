// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef mean_mapping_quality_hpp
#define mean_mapping_quality_hpp

#include <string>
#include <vector>

#include "measure.hpp"

namespace octopus {

class VcfRecord;

namespace csr {

class MeanMappingQuality : public Measure
{
    bool recalculate_;
    std::unique_ptr<Measure> do_clone() const override;
    ResultType do_evaluate(const VcfRecord& call, const FacetMap& facets) const override;
    ResultCardinality do_cardinality() const noexcept override;
    std::string do_name() const override;
    std::vector<std::string> do_requirements() const override;
public:
    MeanMappingQuality(bool recalculate = true);
};

} // namespace csr
} // namespace octopus


#endif