// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef qual_hpp
#define qual_hpp

#include <string>
#include <memory>
#include <utility>

#include "measure.hpp"

namespace octopus {

class VcfRecord;

namespace csr {

class Qual : public Measure
{
public:
    virtual double operator()(const VcfRecord& call) const override;
    virtual std::string name() const override;
};

} // namespace csr
} // namespace octopus


#endif
