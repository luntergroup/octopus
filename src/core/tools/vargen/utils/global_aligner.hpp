// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef global_aligner_hpp
#define global_aligner_hpp

#include <string>

#include "basics/cigar_string.hpp"

namespace octopus { namespace coretools {

struct Model
{
    using ScoreType = short;
    ScoreType match      =  2;
    ScoreType mismatch   = -3;
    ScoreType gap_open   = -8;
    ScoreType gap_extend = -1;
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
