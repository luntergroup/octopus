// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef __Octopus__option_parser__
#define __Octopus__option_parser__

#include <string>
#include <iosfwd>

#include <boost/program_options.hpp>
#include <boost/optional.hpp>

namespace octopus { namespace options {

using OptionMap = boost::program_options::variables_map;

OptionMap parse_options(int argc, const char** argv);

enum class ContigOutputOrder
{
    LexicographicalAscending, LexicographicalDescending,
    ContigSizeAscending, ContigSizeDescending,
    AsInReferenceIndex, AsInReferenceIndexReversed,
    Unspecified
};

struct ContigPloidy
{
    boost::optional<std::string> sample;
    std::string contig;
    unsigned ploidy;
};

enum class RefCallType { Positional, Blocked };

enum class PhasingLevel { Minimal, Conservative, Aggressive };

std::istream& operator>>(std::istream& in, ContigOutputOrder& coo);
std::ostream& operator<<(std::ostream& os, const ContigOutputOrder& coo);
std::istream& operator>>(std::istream& in, ContigPloidy& cp);
std::ostream& operator<<(std::ostream& os, const ContigPloidy& cp);
std::istream& operator>>(std::istream& in, RefCallType& rct);
std::ostream& operator<<(std::ostream& os, const RefCallType& rct);
std::istream& operator>>(std::istream& in, PhasingLevel& pl);
std::ostream& operator<<(std::ostream& os, const PhasingLevel& pl);

} // namespace Options
} // namespace octopus

#endif /* defined(__Octopus__option_parser__) */
