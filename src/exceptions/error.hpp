// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef error_hpp
#define error_hpp

#include <exception>
#include <string>

namespace octopus {

/**
 Error specifies the interface for all octopus exceptions.
 
 std::string is allowed as we won't try to handle any std::bad_alloc.
 */
class Error : public std::exception
{
public:
    virtual ~Error() override = default;
    
    // Who is responsable for this error occurring?
    std::string type() const;
    
    // Where did the error occur (for debugging)? May not be a stacktrace but should give a hint.
    std::string where() const;
    
    // A detailed explanation of *why* the error happened.
    std::string why() const;
    
    // What can be done to resolve the error?
    std::string help() const;
    
    // Just to satisfy std::exception interface, not intended to be called.
    const char* what() const noexcept override;
    
private:
    virtual std::string do_type() const  = 0;
    virtual std::string do_where() const = 0;
    virtual std::string do_why() const   = 0;
    virtual std::string do_help() const  = 0;
    
    mutable std::string what_; // for what()
};

} // namespace octopus

#endif
