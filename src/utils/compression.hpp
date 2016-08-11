// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef compression_hpp
#define compression_hpp

#include <string>

namespace octopus { namespace utils {

std::string compress(const std::string& data);
std::string decompress(const std::string& data);

struct Compress
{
    std::string operator()(const std::string str) const;
};

struct Decompress
{
    std::string operator()(const std::string str) const;
};

} // namespace utils
} // namespace octopus

#endif
