// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "unimplemented_feature_error.hpp"

#include <utility>

namespace octopus {

UnimplementedFeatureError::UnimplementedFeatureError(std::string feature, std::string where)
: feature_ {std::move(feature)}
, where_ {std::move(where)}
{}

std::string UnimplementedFeatureError::do_why() const
{
    return feature_ + " is not currently implemented";
}

std::string UnimplementedFeatureError::do_help() const
{
    return "ensure the specified path is correct and the file is readable";
}

std::string UnimplementedFeatureError::do_where() const
{
    return where_;
}

} // namespace octopus
