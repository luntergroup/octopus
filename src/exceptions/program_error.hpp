// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <string>

#include "error.hpp"

#ifndef program_error_h
#define program_error_h

namespace octopus {

/**
 A ProgramError is any error that octopus is responsible for. Note the error itself may not be a
 bug, but in these cases they should be handled.
 */
class ProgramError : public Error
{
    virtual std::string type() const override { return "program"; }
public:
    virtual ~ProgramError() = default;
};

} // namespace octopus

#endif
