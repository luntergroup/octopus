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
    using Sample = Octopus::SampleIdType;
    
    // Use these to make construction order explicit
    struct Mother { Sample name; };
    struct Father { Sample name; };
    struct Child  { Sample name; };
    
    Trio() = default;
    
    explicit Trio(Mother mother, Father father, Child child);
    
    ~Trio() = default;
    
    Trio(const Trio&)            = default;
    Trio& operator=(const Trio&) = default;
    Trio(Trio&&)                 = default;
    Trio& operator=(Trio&&)      = default;
    
    bool is_child(const Sample& member) const noexcept;
    bool is_parent(const Sample& member) const noexcept;
    
    const Sample& mother() const noexcept;
    const Sample& father() const noexcept;
    const Sample& child() const noexcept;
    
private:
    Sample mother_, father_, child_;
};

#endif /* trio_hpp */
