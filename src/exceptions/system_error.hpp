// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef system_error_h
#define system_error_h

#include <string>

#include "error.hpp"

namespace octopus {

/**
 A SystemError is any error that is not directly attributable to either the user or the program,
 e.g. if a file goes missing or we run out of memory.
 */
class SystemError : public Error
{
    virtual std::string do_type() const override { return "system"; }
public:
    virtual ~SystemError() = default;
};

} // namespace octopus

#endif
