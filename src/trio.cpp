//
//  trio.cpp
//  Octopus
//
//  Created by Daniel Cooke on 29/01/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "trio.hpp"

#include <utility>
#include <stdexcept>

// public member methods

Trio::Trio(Mother mother, Father father, Child child)
:
mother_ {std::move(mother.name)},
father_ {std::move(father.name)},
child_  {std::move(child.name)}
{
    if (mother_ == father_ || mother_ == child_ || father_ == child_) {
        throw std::logic_error {"Trio: members of a trio must be unique"};
    }
}

bool Trio::is_child(const Sample& member) const noexcept
{
    return member == child_;
}

bool Trio::is_parent(const Sample& member) const noexcept
{
    return member == mother_ || member == father_;
}

const Trio::Sample& Trio::mother() const noexcept
{
    return mother_;
}

const Trio::Sample& Trio::father() const noexcept
{
    return father_;
}

const Trio::Sample& Trio::child() const noexcept
{
    return child_;
}
