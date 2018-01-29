// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "unwritable_file_error.hpp"

#include <utility>
#include <sstream>

namespace octopus {

UnwritableFileError::UnwritableFileError(Path file) : file_ {std::move(file)} {}

UnwritableFileError::UnwritableFileError(Path file, std::string type)
: file_ {std::move(file)}
, type_ {std::move(type)}
{}

std::string UnwritableFileError::do_why() const
{
    std::ostringstream ss {};
    ss << "the ";
    if (type_) {
        ss << *type_ << ' ';
    }
    ss << "file you specified " << file_ << ' ';
    ss << "is not writable";
    return ss.str();
}

std::string UnwritableFileError::do_help() const
{
    return "ensure the specified path is correct and the location is writable (check permissions)";
}

} // namespace octopus
