// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef option_parser_hpp
#define option_parser_hpp

#include <string>
#include <iosfwd>

#include <boost/program_options.hpp>
#include <boost/optional.hpp>

namespace octopus { namespace options {

using OptionMap = boost::program_options::variables_map;

OptionMap parse_options(int argc, const char** argv);

enum class ContigOutputOrder
{
    lexicographicalAscending, lexicographicalDescending,
    contigSizeAscending, contigSizeDescending,
    asInReferenceIndex, asInReferenceIndexReversed,
    unspecified
};

struct ContigPloidy
{
    boost::optional<std::string> sample;
    std::string contig;
    int ploidy;
};

enum class RefCallType { positional, blocked };

enum class PhasingLevel { minimal, conservative, aggressive };

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

#endif
