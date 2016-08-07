// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "missing_file_error.hpp"

#include <utility>
#include <sstream>

namespace octopus {

MissingFileError::MissingFileError(Path file, std::string type)
: file_ {std::move(file)}, type_ {std::move(type)} {}

std::string MissingFileError::do_why() const
{
    std::ostringstream ss {};
    
    ss << "the " << type_ << " file you specified " << file_ << " does not exist";
    
    return ss.str();
}

std::string MissingFileError::do_help() const
{
    return "ensure the specified path is correct and the file is readable";
}

} // namespace octopus
