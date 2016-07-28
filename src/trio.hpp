//
//  trio.hpp
//  Octopus
//
//  Created by Daniel Cooke on 29/01/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef trio_hpp
#define trio_hpp

#include "common.hpp"

class Trio
{
public:
    using SampleName = octopus::SampleName;
    
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

#endif /* trio_hpp */
