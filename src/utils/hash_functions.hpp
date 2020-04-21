// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef hash_functions_hpp
#define hash_functions_hpp

#include <string>

#include <boost/utility/string_ref.hpp>
#include <boost/functional/hash.hpp>
#include <boost/filesystem/path.hpp>

namespace octopus { namespace utils
{
    struct StringRefHash
    {
        std::size_t operator()(const boost::string_ref& str) const
        {
            return boost::hash_range(str.begin(), str.end());
        }
    };
    
    struct FilepathHash
    {
        std::size_t operator()(const boost::filesystem::path& path) const
        {
            return std::hash<std::string>()(path.string());
        }
    };
} // namespace utils
} // namespace octopus

#endif
