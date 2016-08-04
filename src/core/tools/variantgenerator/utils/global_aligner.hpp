//
//  global_aligner.hpp
//  octopus
//
//  Created by Daniel Cooke on 03/07/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef global_aligner_hpp
#define global_aligner_hpp

#include <string>

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
    std::string cigar;
    int score;
};

Alignment align(const std::string& target, const std::string& query, Model model = Model {});

} // namespace coretools
} // namespace octopus


#endif /* global_aligner_hpp */
