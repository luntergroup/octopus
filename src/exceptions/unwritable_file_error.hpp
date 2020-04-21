// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef unwritable_file_error_hpp
#define unwritable_file_error_hpp

#include <string>

#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>

#include "user_error.hpp"

namespace octopus {

class UnwritableFileError : public UserError
{
public:
    using Path = boost::filesystem::path;
    
    UnwritableFileError() = delete;
    
    UnwritableFileError(Path file);
    
    UnwritableFileError(Path file, std::string type);
    
    virtual ~UnwritableFileError() override = default;
    
private:
    virtual std::string do_why() const override;
    virtual std::string do_help() const override;
    
    Path file_;
    boost::optional<std::string> type_;
};

} // namespace octopus

#endif
