//
//  compression.hpp
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_compression_hpp
#define Octopus_compression_hpp

#include <string>

namespace Octopus
{
    std::string compress(const std::string& data);
    std::string decompress(const std::string& data);
    
    struct Compress
    {
        std::string operator()(const std::string str) const;
    };
    
    struct Decompress
    {
        std::string operator()(const std::string str) const;
    };
} // namespace Octopus

#endif
