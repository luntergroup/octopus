// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef config_hpp
#define config_hpp

#include <string>
#include <vector>
#include <iosfwd>

#include <boost/optional.hpp>

namespace octopus { namespace config {

struct VersionNumber
{
    unsigned short major, minor;
    boost::optional<unsigned short> patch = boost::none;
    boost::optional<std::string> name = boost::none;
};

extern const VersionNumber Version;

std::ostream& operator<<(std::ostream& os, const VersionNumber& version);

extern const std::string HelpForum, BugReport;

extern const std::vector<std::string> Authors;

extern const std::string CopyrightNotice;

extern const unsigned CommandLineWidth;

} // namespace config
} // namespace octopus

#endif
