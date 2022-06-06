// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "config.hpp"

#include <ostream>
#include <sstream>
#include <algorithm>

#include <boost/range/algorithm_ext/erase.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "version.hpp"
#include "system.hpp"

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

static auto get_simd_extension()
{
    #if defined(__AVX2__)
        return SystemInfo::SIMDExtension::avx512;
    #elif defined(__AVX512F__) && defined(__AVX512BW__)
        return SystemInfo::SIMDExtension::avx2;
    #elif defined(__SSE2__)
        return SystemInfo::SIMDExtension::sse2;
    #else
        return SystemInfo::SIMDExtension::neon;
    #endif
}

const SystemInfo System {SYSTEM_PROCESSOR,
                         SYSTEM_NAME,
                         SYSTEM_VERSION,
                         COMPILER_NAME,
                         COMPILER_VERSION,
                         BOOSTLIB_VERSION,
                         BUILD_TYPE,
                         get_simd_extension()};

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

std::ostream& operator<<(std::ostream& os, const SystemInfo::SIMDExtension simd)
{
    using SIMD = SystemInfo::SIMDExtension;
    switch (simd) {
        case SIMD::sse2: os << "SSE2"; break;
        case SIMD::avx2: os << "AVX2"; break;
        case SIMD::avx512: os << "AVX512"; break;
        case SIMD::neon: os << "NEON"; break;
    }
    return os;
}

std::string to_string(const VersionNumber& version, const bool spaces)
{
    std::ostringstream ss {};
    ss << config::Version;
    auto result = ss.str();
    if (!spaces) {
        std::replace(result.begin(), result.end(), ' ', '_');
        boost::remove_erase_if(result, boost::is_any_of("()"));
    }
    return result;
}

const std::string HelpForum {"https://github.com/luntergroup/octopus/issues"};

const std::string BugReport {"https://github.com/luntergroup/octopus/issues"};

const std::vector<std::string> Authors {"Daniel Cooke"};

const std::string CopyrightNotice {"Copyright (c) 2015-2021 University of Oxford"};

const unsigned CommandLineWidth {72};

} // namespace config
} // namespace octopus
