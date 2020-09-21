// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "path_utils.hpp"

#include <string>
#include <sstream>

#include <boost/filesystem/operations.hpp>

#include "exceptions/system_error.hpp"

namespace octopus {

boost::optional<fs::path> get_home_directory()
{
    static const auto env = std::getenv("HOME");
    if (env == nullptr) return boost::none;
    const fs::path home {env};
    if (fs::is_directory(home)) return home;
    return boost::none;
}

bool is_shorthand_user_path(const fs::path& path) noexcept
{
    return !path.empty() && path.string().front() == '~';
}

class UnknownHomeDirectory : public SystemError
{
    std::string do_where() const override
    {
        return "expand_user_path";
    }

    std::string do_why() const override
    {
        std::ostringstream ss {};
        ss << "Unable to expand shorthand path you specified ";
        ss << path_;
        ss << " as your home directory cannot be located";
        return ss.str();
    }

    std::string do_help() const override
    {
        return "ensure your HOME environment variable is set properly";
    }

    fs::path path_;
public:
    UnknownHomeDirectory(fs::path p) : path_ {std::move(p)} {}
};

fs::path expand_user_path(const fs::path& path)
{
    if (is_shorthand_user_path(path)) {
        if (path.string().size() > 1 && path.string()[1] == '/') {
            const auto home_dir = get_home_directory();
            if (home_dir) {
                return fs::path {home_dir->string() + path.string().substr(1)};
            }
            throw UnknownHomeDirectory {path};
        }
        return path;
    }
    return path;
}

fs::path
resolve_path(const fs::path& path, const fs::path& working_directory,
             const WorkingDirectoryResolvePolicy wd_policy,
             const SymblinkResolvePolicy symlink_policy)
{
    fs::path result {};
    if (is_shorthand_user_path(path)) {
        result = expand_user_path(path); // must be a root path
    } else {
        const auto wd_full_path = working_directory / path;
        if (fs::exists(path)) {
            if (fs::exists(wd_full_path) && wd_policy == WorkingDirectoryResolvePolicy::prefer_working_directory) {
                result = wd_full_path;
            } else {
                result = path; // must be a root path
            }
        } else {
            const auto parent_dir = path.parent_path();
            if (fs::exists(parent_dir) && fs::is_directory(parent_dir)) {
                auto tmp = working_directory;
                tmp /= path;
                auto wd_parent = tmp.parent_path();
                if (fs::exists(wd_parent) && fs::is_directory(wd_parent)) {
                    result = tmp; // prefer working directory in case of name clash
                } else {
                    result = path; // must be yet-to-be-created root path
                }
            } else {
                result = wd_full_path;
            }
        }
    }
    if (symlink_policy == SymblinkResolvePolicy::resolve) {
        result = fs::canonical(result);
    } else {
        result = fs::absolute(result);
    }
    return result;
}

} // namespace octopus
