// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef missing_index_error_hpp
#define missing_index_error_hpp

#include <string>

#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>

#include "user_error.hpp"

namespace octopus {

/**
 A MissingIndexError should be thrown when an index associated with a user specified file
 does not exist.
 */
class MissingIndexError : public UserError
{
public:
    using Path = boost::filesystem::path;
    
    MissingIndexError() = delete;
    
    // When no associated index could be located
    MissingIndexError(Path associate, std::string type);
    
    // When an index was given, but is non-existent
    MissingIndexError(Path associate, Path given_index, std::string type);
    
    virtual ~MissingIndexError() override = default;
    
private:
    virtual std::string do_why() const override;
    virtual std::string do_help() const override;
    
    Path associate_;
    boost::optional<Path> given_index_;
    std::string type_;
};

} // namespace octopus

#endif
