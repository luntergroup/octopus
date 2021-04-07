// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef missing_file_error_hpp
#define missing_file_error_hpp

#include <string>

#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>

#include "user_error.hpp"

namespace octopus {

/**
 A MissingFileError should be thrown when a user-specified file does not exist.
 */
class MissingFileError : public UserError
{
public:
    using Path = boost::filesystem::path;
    
    MissingFileError() = delete;
    
    MissingFileError(Path file);
    
    MissingFileError(Path file, std::string type);
    
    virtual ~MissingFileError() override = default;
    
    void set_location_specified(std::string location) noexcept;
    
private:
    virtual std::string do_why() const override;
    virtual std::string do_help() const override;
    
    Path file_;
    boost::optional<std::string> type_, location_;
};

} // namespace octopus

#endif
