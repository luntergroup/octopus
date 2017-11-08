// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef trio_hpp
#define trio_hpp

#include "config/common.hpp"

namespace octopus {

class Trio
{
public:
    // Use these to make construction order explicit
    struct Mother { SampleName name; };
    struct Father { SampleName name; };
    struct Child  { SampleName name; };
    
    Trio() = default;
    
    Trio(Mother mother, Father father, Child child);
    
    Trio(const Trio&)            = default;
    Trio& operator=(const Trio&) = default;
    Trio(Trio&&)                 = default;
    Trio& operator=(Trio&&)      = default;
    
    ~Trio() = default;
    
    bool is_child(const SampleName& member) const noexcept;
    bool is_parent(const SampleName& member) const noexcept;
    
    const SampleName& mother() const noexcept;
    const SampleName& father() const noexcept;
    const SampleName& child() const noexcept;
    
private:
    SampleName mother_, father_, child_;
};

} // namespace octopus

#endif
