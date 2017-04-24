// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "config.hpp"

#include <ostream>

namespace octopus { namespace config {

const VersionNumber Version {0, 2, 0, boost::optional<std::string> {"alpha"}};

std::ostream& operator<<(std::ostream& os, const VersionNumber& version)
{
    os << version.major << '.' << version.minor;
    if (version.patch) os << '.' << *version.patch;
    if (version.name) os << '-' << *version.name;
    return os;
}

const std::string HelpForum {"https://github.com/luntergroup/octopus/issues"};

const std::string BugReport {"https://github.com/luntergroup/octopus/issues"};

const std::vector<std::string> Authors {"Daniel Cooke"};

const std::string CopyrightNotice {"Copyright (c) 2017 University of Oxford"};

const unsigned CommandLineWidth {72};

} // namespace config
} // namespace octopus
