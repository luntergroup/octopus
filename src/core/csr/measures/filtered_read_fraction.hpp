// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef filtered_read_fraction_hpp
#define filtered_read_fraction_hpp

#include <string>
#include <vector>

#include "measure.hpp"
#include "depth.hpp"

namespace octopus {

class VcfRecord;

namespace csr {

class FilteredReadFraction : public Measure
{
    Depth calling_depth_, filtering_depth_;
    std::unique_ptr<Measure> do_clone() const override;
    ResultType do_evaluate(const VcfRecord& call, const FacetMap& facets) const override;
    std::string do_name() const override;
    std::vector<std::string> do_requirements() const override;
public:
    FilteredReadFraction();
};

} // namespace csr
} // namespace octopus

#endif
