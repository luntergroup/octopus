// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "file_open_error.hpp"

#include <utility>
#include <sstream>

namespace octopus {

FileOpenError::FileOpenError(Path file) : file_ {std::move(file)} {}

FileOpenError::FileOpenError(Path file, std::string type)
: file_ {std::move(file)}
, type_ {std::move(type)}
{}

std::string FileOpenError::do_why() const
{
    std::ostringstream ss {};
    ss << "the ";
    if (type_) {
        ss << *type_ << ' ';
    }
    ss << "file " << file_ << " could not be opened. The OS open file limit (ulimit -n) may have been exceeded";
    return ss.str();
}

std::string FileOpenError::do_help() const
{
    return "Increase the OS open file limit";
}

std::string FileOpenError::do_where() const
{
    return "fopen";
}
    
} // namespace octopus
