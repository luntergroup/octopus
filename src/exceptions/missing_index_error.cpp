// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "missing_index_error.hpp"

#include <utility>
#include <sstream>

namespace octopus {

MissingIndexError::MissingIndexError(Path associate, std::string type)
: associate_ {std::move(associate)}
, type_ {std::move(type)}
{}

MissingIndexError::MissingIndexError(Path associate, Path given_index, std::string type)
: associate_ {std::move(associate)}
, given_index_ {std::move(given_index)}
, type_ {std::move(type)}
{}

std::string MissingIndexError::do_why() const
{
    std::ostringstream ss {};
    if (given_index_) {
        ss << "The index file that you specified " << *given_index_ << " for the " << type_ << " file "
            << associate_ << " does not exist";
    } else {
        ss << "No associated index file could be found for the " << type_ << " file " << associate_;
    }
    return ss.str();
}

std::string MissingIndexError::do_help() const
{
    return "ensure the specified path is correct and the file is readable";
}

} // namespace octopus
