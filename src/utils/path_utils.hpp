// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef path_utils_hpp
#define path_utils_hpp

#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>

namespace octopus {

namespace fs = boost::filesystem;

boost::optional<fs::path> get_home_directory();

bool is_shorthand_user_path(const fs::path& path) noexcept;

fs::path expand_user_path(const fs::path& path);

fs::path resolve_path(const fs::path& path, const fs::path& working_directory);

} // namespace octopus

#endif
