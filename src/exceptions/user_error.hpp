// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef user_error_hpp
#define user_error_hpp

#include <string>

#include "error.hpp"

namespace octopus {

/**
 A UserError is used for any error caused by bad user input, e.g. missing files, incorrect file types,
 or invalid command line values.
 */
class UserError : public Error
{
    virtual std::string do_type() const override { return "user"; }
public:
    virtual ~UserError() = default;
};

} // namespace octopus

#endif
