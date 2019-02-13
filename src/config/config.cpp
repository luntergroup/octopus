// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "config.hpp"

#include <ostream>

#include "version.hpp"

namespace octopus { namespace config {

static boost::optional<std::string> get_release_name()
{
    std::string name {VERSION_RELEASE};
    if (name.empty()) {
        return boost::none;
    } else {
        return name;
    }
}

static boost::optional<std::string> get_git_branch_name()
{
    std::string name {GIT_BRANCH};
    if (name.empty()) {
        return boost::none;
    } else {
        return name;
    }
}

static boost::optional<std::string> get_git_commit()
{
    std::string name {GIT_COMMIT_HASH};
    if (name.empty()) {
        return boost::none;
    } else {
        return name;
    }
}

const VersionNumber Version {VERSION_MAJOR,
                             VERSION_MINOR,
                             VERSION_PATCH,
                             get_release_name(),
                             get_git_branch_name(),
                             get_git_commit()};

std::ostream& operator<<(std::ostream& os, const VersionNumber& version)
{
    os << version.major << '.' << version.minor;
    if (version.patch) os << '.' << *version.patch;
    if (version.name) os << '-' << *version.name;
    if ((version.branch && *version.branch != "master") || version.commit) {
        os << " (";
        std::string space {""};
        if (version.branch && *version.branch != "master") {
            os << *version.branch;
            space = " ";
        }
        if (version.commit) os << space << *version.commit;
        os << ')';
    }
    return os;
}

const std::string HelpForum {"https://github.com/luntergroup/octopus/issues"};

const std::string BugReport {"https://github.com/luntergroup/octopus/issues"};

const std::vector<std::string> Authors {"Daniel Cooke"};

const std::string CopyrightNotice {"Copyright (c) 2015-2019 University of Oxford"};

const unsigned CommandLineWidth {72};

} // namespace config
} // namespace octopus
