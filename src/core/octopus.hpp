// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef octopus_hpp
#define octopus_hpp

#include <string>

#include "calling_components.hpp"

namespace octopus {

struct UserCommandInfo
{
    std::string command, options;
};

void run_octopus(GenomeCallingComponents& components, UserCommandInfo info);

}

#endif
