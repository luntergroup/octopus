// Copyright (c) 2015-2021 Daniel Cooke
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

FileOpenError::FileOpenError(Path file, std::error_code error)
: file_ {std::move(file)}
, error_ {std::move(error)}
{}

FileOpenError::FileOpenError(Path file, std::string type, std::error_code error)
: file_ {std::move(file)}
, type_ {std::move(type)}
, error_ {std::move(error)}
{}

std::string FileOpenError::do_why() const
{
    std::ostringstream ss {};
    ss << "the ";
    if (type_) {
        ss << *type_ << ' ';
    }
    ss << "file " << file_ << " could not be opened";
    if (error_) {
        ss << ": " << error_->message();
    }
    return ss.str();
}

std::string FileOpenError::do_help() const
{
    return "Check the output path is writable";
}

std::string FileOpenError::do_where() const
{
    return "fopen";
}
    
} // namespace octopus
