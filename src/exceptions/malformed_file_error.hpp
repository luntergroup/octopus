// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef malformed_file_error_hpp
#define malformed_file_error_hpp

#include <string>
#include <vector>

#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>

#include "user_error.hpp"

namespace octopus {

/**
 A MalformedFileError should be thrown when a user-specified file is of the wrong type,
 or is corrupted.
 */
class MalformedFileError : public UserError
{
public:
    using Path = boost::filesystem::path;
    
    MalformedFileError() = delete;
    
    MalformedFileError(Path file);
    
    MalformedFileError(Path file, std::string required_type);
    
    MalformedFileError(Path file, std::vector<std::string> valid_types);
    
    virtual ~MalformedFileError() override = default;
    
    void set_reason(std::string reason) noexcept;
    
    void set_location_specified(std::string location) noexcept;
    
private:
    virtual std::string do_why() const override;
    virtual std::string do_help() const override;
    
    Path file_;
    std::vector<std::string> valid_types_;
    boost::optional<std::string> reason_, location_;
};

} // namespace octopus

#endif
