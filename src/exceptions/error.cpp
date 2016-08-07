// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "error.hpp"

#include <sstream>

namespace octopus {

std::string Error::type() const
{
    return do_type();
}

std::string Error::where() const
{
    return do_where();
}

std::string Error::why() const
{
    return do_why();
}

std::string Error::help() const
{
    return do_help();
}

const char* Error::what() const noexcept
{
    std::ostringstream ss {};
    
    what_ = ss.str();
    
    ss << "type: " << type();
    
    return what_.c_str();
}

} // namespace octopus
