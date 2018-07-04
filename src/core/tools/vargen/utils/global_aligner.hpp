// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef global_aligner_hpp
#define global_aligner_hpp

#include <string>

#include "basics/cigar_string.hpp"

namespace octopus { namespace coretools {

struct Model
{
    short match      =  2;
    short mismatch   = -3;
    short gap_open   = -8;
    short gap_extend = -1;
};

struct Alignment
{
    CigarString cigar;
    int score;
};

Alignment align(const std::string& target, const std::string& query, Model model = Model {});

} // namespace coretools
} // namespace octopus


#endif
