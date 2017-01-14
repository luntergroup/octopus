// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef measure_hpp
#define measure_hpp

#include <string>
#include <unordered_set>
#include <memory>
#include <utility>

//#include "../facets/facet.hpp"

namespace octopus {

class VcfRecord;

namespace csr {

class Measure
{
public:
    virtual double operator()(const VcfRecord& call) const = 0;
    virtual std::string name() const = 0;
    //virtual std::unordered_set<Facet> requirements() const noexcept { return {}; }
};

class MeasureWrapper
{
public:
    MeasureWrapper() = delete;
    
    MeasureWrapper(std::unique_ptr<Measure> measure) : measure_ {std::move(measure)} {}
    
    MeasureWrapper(const MeasureWrapper&)            = delete;
    MeasureWrapper& operator=(const MeasureWrapper&) = delete;
    MeasureWrapper(MeasureWrapper&&)                 = default;
    MeasureWrapper& operator=(MeasureWrapper&&)      = default;
    
    ~MeasureWrapper() = default;
    
    auto operator()(const VcfRecord& call) const { return (*measure_)(call); }
    std::string name() const { return measure_->name(); }
    
private:
    std::unique_ptr<Measure> measure_;
};

template <typename M, typename... Args>
MeasureWrapper make_wrapped_measure(Args&&... args)
{
    return MeasureWrapper {std::make_unique<M>(std::forward<Args>(args)...)};
}

} // namespace csr
} // namespace octopus

#endif
