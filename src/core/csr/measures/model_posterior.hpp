// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef model_posterior_hpp
#define model_posterior_hpp

#include <string>

#include "measure.hpp"

namespace octopus {

class VcfRecord;

namespace csr {

class ModelPosterior : public Measure
{
    std::unique_ptr<Measure> do_clone() const override;
    ResultType do_evaluate(const VcfRecord& call, const FacetMap& facets) const override;
    std::string do_name() const override;
};

} // namespace csr
} // namespace octopus


#endif
