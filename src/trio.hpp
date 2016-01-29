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
    using Member = Octopus::SampleIdType;
    
    // Use these to make construction order explicit
    struct Mother { Member name; };
    struct Father { Member name; };
    struct Child  { Member name; };
    
    Trio() = default;
    explicit Trio(Mother mother, Father father, Child child);
    ~Trio() = default;
    
    Trio(const Trio&)            = default;
    Trio& operator=(const Trio&) = default;
    Trio(Trio&&)                 = default;
    Trio& operator=(Trio&&)      = default;
    
    bool is_child(const Member& member) const noexcept;
    bool is_parent(const Member& member) const noexcept;
    
    const Member& get_mother() const noexcept;
    const Member& get_father() const noexcept;
    const Member& get_child() const noexcept;
    
private:
    Member mother_, father_, child_;
};

#endif /* trio_hpp */
