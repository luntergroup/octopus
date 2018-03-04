// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <string>

#include <boost/filesystem/path.hpp>

#include "user_error.hpp"

#ifndef unwritable_file_error_hpp
#define unwritable_file_error_hpp

namespace octopus {

class UnwritableFileError : public UserError
{
    
};

} // namespace octopus

#endif
