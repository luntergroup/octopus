// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <string>

#include "error.hpp"
#include "config/config.hpp"

#ifndef program_error_hpp
#define program_error_hpp

namespace octopus {

/**
 A ProgramError is any error that octopus is responsible for. Note the error itself may not be a
 bug, but in these cases they should be handled.
 */
class ProgramError : public Error
{
    virtual std::string do_type() const override { return "program"; }
    virtual std::string do_help() const override
    {
        return "Run in debug mode and send the log file to " + config::BugReport;
    }
public:
    virtual ~ProgramError() = default;
};

} // namespace octopus

#endif
