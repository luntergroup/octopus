//
//  option_parser.hpp
//  Octopus
//
//  Created by Daniel Cooke on 27/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__option_parser__
#define __Octopus__option_parser__

#include <string>
#include <iosfwd>

#include <boost/program_options.hpp>
#include <boost/optional.hpp>

namespace Octopus
{
namespace Options
{
    using OptionMap = boost::program_options::variables_map;
    
    boost::optional<OptionMap> parse_options(int argc, const char** argv);
    
    bool is_run_command(const OptionMap& options);
    
    enum class ContigOutputOrder
    {
        LexicographicalAscending, LexicographicalDescending,
        ContigSizeAscending, ContigSizeDescending,
        AsInReferenceIndex, AsInReferenceIndexReversed,
        Unspecified
    };
    
    struct ContigPloidy
    {
        std::string contig;
        unsigned ploidy;
    };
    
    enum class RefCallType { Positional, Blocked };
    
    enum class PhasingLevel { Minimal, Conservative, Aggressive };
    
    std::istream& operator>>(std::istream& in, ContigPloidy& cp);
    std::ostream& operator<<(std::ostream& os, const ContigPloidy& cp);
} // namespace Options
} // namespace Octopus

#endif /* defined(__Octopus__option_parser__) */
