// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef malformed_file_error_hpp
#define malformed_file_error_hpp

#include <string>
#include <vector>

#include <boost/filesystem/path.hpp>

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
    
    MalformedFileError(Path name);
    
    MalformedFileError(Path name, std::string required_type);
    
    MalformedFileError(Path name, std::vector<std::string> valid_types);
    
    virtual ~MalformedFileError() override = default;
    
private:
    virtual std::string do_why() const override;
    virtual std::string do_help() const override;
    
    Path name_;
    std::vector<std::string> valid_types_;
};

} // namespace octopus

#endif
