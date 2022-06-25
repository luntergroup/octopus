// Copyright (c) 2015-2021 Daniel Cooke
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
    referenceIndex, referenceIndexReversed,
    unspecified
};

struct ContigPloidy
{
    boost::optional<std::string> sample;
    std::string contig;
    int ploidy;
};

enum class ExtensionLevel { minimal, conservative, moderate, aggressive, unlimited };
enum class LaggingLevel { none, conservative, moderate, optimistic, aggressive };
enum class BacktrackLevel { none, moderate, aggressive };
enum class NormalContaminationRisk { low, high };
enum class BadRegionTolerance { low, normal, high, unlimited };
enum class ReadLinkage { none, paired, linked };
enum class CandidateVariantDiscoveryProtocol { illumina, pacbio };
enum class RealignedBAMType { full, mini };
enum class ReadDeduplicationDetectionPolicy { relaxed, aggressive };
enum class ModelPosteriorPolicy { all, off, special };
enum class PhasingPolicy { conservative, aggressive, automatic };

struct SampleDropoutConcentrationPair
{
    std::string sample;
    float concentration;
};

struct SamTag
{
    std::string tag;
    boost::optional<std::string> value;
};

std::istream& operator>>(std::istream& in, ContigOutputOrder& order);
std::ostream& operator<<(std::ostream& os, const ContigOutputOrder& order);
std::istream& operator>>(std::istream& in, ContigPloidy& plodies);
std::ostream& operator<<(std::ostream& os, const ContigPloidy& plodies);
std::istream& operator>>(std::istream& in, ExtensionLevel& level);
std::ostream& operator<<(std::ostream& os, const ExtensionLevel& level);
std::istream& operator>>(std::istream& in, BacktrackLevel& level);
std::ostream& operator<<(std::ostream& os, const BacktrackLevel& level);
std::istream& operator>>(std::istream& in, LaggingLevel& level);
std::ostream& operator<<(std::ostream& os, const LaggingLevel& level);
std::istream& operator>>(std::istream& in, NormalContaminationRisk& risk);
std::ostream& operator<<(std::ostream& os, const NormalContaminationRisk& risk);
std::istream& operator>>(std::istream& in, BadRegionTolerance& risk);
std::ostream& operator<<(std::ostream& os, const BadRegionTolerance& risk);
std::istream& operator>>(std::istream& in, ReadLinkage& linkage);
std::ostream& operator<<(std::ostream& os, const ReadLinkage& linkage);
std::istream& operator>>(std::istream& in, CandidateVariantDiscoveryProtocol& protocol);
std::ostream& operator<<(std::ostream& os, const CandidateVariantDiscoveryProtocol& protocol);
std::istream& operator>>(std::istream& in, RealignedBAMType& type);
std::ostream& operator<<(std::ostream& os, const RealignedBAMType& type);
std::istream& operator>>(std::istream& in, ReadDeduplicationDetectionPolicy& type);
std::ostream& operator<<(std::ostream& os, const ReadDeduplicationDetectionPolicy& type);
std::istream& operator>>(std::istream& in, ModelPosteriorPolicy& policy);
std::ostream& operator<<(std::ostream& os, const ModelPosteriorPolicy& policy);
std::istream& operator>>(std::istream& in, SampleDropoutConcentrationPair& concentration);
std::ostream& operator<<(std::ostream& os, const SampleDropoutConcentrationPair& concentration);
std::istream& operator>>(std::istream& in, PhasingPolicy& policy);
std::ostream& operator<<(std::ostream& os, const PhasingPolicy& policy);
std::istream& operator>>(std::istream& in, SamTag& policy);
std::ostream& operator<<(std::ostream& os, const SamTag& policy);

std::ostream& operator<<(std::ostream& os, const OptionMap& options);
std::string to_string(const OptionMap& options, bool one_line = false, bool mark_modified = true);

} // namespace options
} // namespace octopus

#endif
