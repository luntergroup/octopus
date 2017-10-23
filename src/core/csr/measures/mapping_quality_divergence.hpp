// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef mapping_quality_divergence_hpp
#define mapping_quality_divergence_hpp

#include <string>

#include "measure.hpp"

namespace octopus {

class VcfRecord;

namespace csr {

class MappingQualityDivergence : public Measure
{
    ResultType do_evaluate(const VcfRecord& call, const FacetMap& facets) const override;
    std::string do_name() const override;
    std::vector<std::string> do_requirements() const override;
};

} // namespace csr
} // namespace octopus

#endif