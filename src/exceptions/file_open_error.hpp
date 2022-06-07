// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef file_open_error_hpp
#define file_open_error_hpp

#include <string>

#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>

#include "program_error.hpp"

namespace octopus {

class FileOpenError : public ProgramError
{
public:
    using Path = boost::filesystem::path;
    
    FileOpenError() = delete;
    
    FileOpenError(Path file);
    FileOpenError(Path file, std::string type);
    FileOpenError(Path file, std::error_code error);
    FileOpenError(Path file, std::string type, std::error_code error);
    
    virtual ~FileOpenError() override = default;

private:
    virtual std::string do_why() const override;
    virtual std::string do_help() const override;
    virtual std::string do_where() const override;
    
    Path file_;
    boost::optional<std::string> type_;
    boost::optional<std::error_code> error_;
};
    
} // namespace octopus

#endif
