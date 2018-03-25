// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef quality_hpp
#define quality_hpp

#include <string>

#include "measure.hpp"

namespace octopus {

class VcfRecord;

namespace csr {

class Quality : public Measure
{
    std::unique_ptr<Measure> do_clone() const override;
    ResultType do_evaluate(const VcfRecord& call, const FacetMap& facets) const override;
    ResultCardinality do_cardinality() const noexcept override;
    std::string do_name() const override;
    std::string do_describe() const override;
    bool is_required_vcf_field() const noexcept override { return true; }
};

} // namespace csr
} // namespace octopus

#endif
