// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "trio.hpp"

#include <utility>
#include <stdexcept>

namespace octopus {

// public member methods

Trio::Trio(Mother mother, Father father, Child child)
: mother_ {std::move(mother.name)}
, father_ {std::move(father.name)}
, child_  {std::move(child.name)}
{
    if (mother_ == father_ || mother_ == child_ || father_ == child_) {
        throw std::logic_error {"Trio: members of a trio must be unique"};
    }
}

bool Trio::is_child(const SampleName& member) const noexcept
{
    return member == child_;
}

bool Trio::is_parent(const SampleName& member) const noexcept
{
    return member == mother_ || member == father_;
}

const SampleName& Trio::mother() const noexcept
{
    return mother_;
}

const SampleName& Trio::father() const noexcept
{
    return father_;
}

const SampleName& Trio::child() const noexcept
{
    return child_;
}

} // namespace octopus
