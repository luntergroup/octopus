// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef error_handler_hpp
#define error_handler_hpp

#include <exception>
#include <new>

#include "exceptions/error.hpp"

namespace octopus {

void log_error(const Error& error);
void log_error(const std::bad_alloc& error);
void log_error(const std::exception& error);

void log_unknown_error();

} // namepace octopus

#endif
