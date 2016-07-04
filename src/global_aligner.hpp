//
//  global_aligner.hpp
//  Octopus
//
//  Created by Daniel Cooke on 03/07/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef global_aligner_hpp
#define global_aligner_hpp

#include <string>
#include <tuple>

struct Model
{
    short match      =  2;
    short mismatch   = -3;
    short gap_open   = -8;
    short gap_extend = -1;
};

std::pair<std::string, int>
align(const std::string& target, const std::string& query, Model model = Model {});

#endif /* global_aligner_hpp */
