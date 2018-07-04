// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "system_utils.hpp"

#include <sys/resource.h>

namespace octopus {

std::size_t get_max_open_files()
{
    struct rlimit lim;
    getrlimit(RLIMIT_NOFILE, &lim);
    return lim.rlim_cur;
}

} // namespace octopus
