// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef quality_by_depth_hpp
#define quality_by_depth_hpp

#include <string>
#include <memory>
#include <utility>

#include "measure.hpp"

namespace octopus {

class VcfRecord;

namespace csr
{
    class QualityByDepth : public Measure
    {
        virtual double operator()(const VcfRecord& call) const override;
        virtual std::string name() const override;
    };
} // namespace csr
} // namespace octopus

#endif /* quality_by_depth_hpp */
