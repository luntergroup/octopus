// Copyright (c) 2015-2019 Daniel Cooke
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
enum class ExtensionLevel { conservative, normal, optimistic, aggressive };
enum class LaggingLevel { minimal, conservative, moderate, normal, aggressive };
enum class BacktrackLevel { none, normal, aggressive };
enum class NormalContaminationRisk { low, high };
enum class CandidateVariantDiscoveryProtocol { illumina, pacbio };

std::istream& operator>>(std::istream& in, ContigOutputOrder& order);
std::ostream& operator<<(std::ostream& os, const ContigOutputOrder& order);
std::istream& operator>>(std::istream& in, ContigPloidy& plodies);
std::ostream& operator<<(std::ostream& os, const ContigPloidy& plodies);
std::istream& operator>>(std::istream& in, RefCallType& type);
std::ostream& operator<<(std::ostream& os, const RefCallType& type);
std::istream& operator>>(std::istream& in, ExtensionLevel& level);
std::ostream& operator<<(std::ostream& os, const ExtensionLevel& level);
std::istream& operator>>(std::istream& in, BacktrackLevel& level);
std::ostream& operator<<(std::ostream& os, const BacktrackLevel& level);
std::istream& operator>>(std::istream& in, LaggingLevel& level);
std::ostream& operator<<(std::ostream& os, const LaggingLevel& level);
std::istream& operator>>(std::istream& in, NormalContaminationRisk& risk);
std::ostream& operator<<(std::ostream& os, const NormalContaminationRisk& risk);
std::istream& operator>>(std::istream& in, CandidateVariantDiscoveryProtocol& protocol);
std::ostream& operator<<(std::ostream& os, const CandidateVariantDiscoveryProtocol& protocol);

} // namespace Options
} // namespace octopus

#endif
