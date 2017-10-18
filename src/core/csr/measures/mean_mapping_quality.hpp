// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef mean_mapping_quality_hpp
#define mean_mapping_quality_hpp

#include <string>
#include <memory>
#include <utility>

#include "measure.hpp"

namespace octopus {

class VcfRecord;

namespace csr {

class MeanMappingQuality : public Measure
{
    virtual ResultType do_evaluate(const VcfRecord& call) const override;
    virtual std::string do_name() const override;
};

} // namespace csr
} // namespace octopus


#endif