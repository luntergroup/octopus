//
//  read_indel_error_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 17/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef read_indel_error_model_hpp
#define read_indel_error_model_hpp

#include <vector>
#include <cstdint>

#include "haplotype.hpp"

namespace Octopus
{
    class ReadIndelErrorModel
    {
    public:
        std::vector<std::uint8_t> calculate_gap_open_penalties(const Haplotype& haplotype) const;
        
    private:
        
    };
} // namespace Octopus

#endif /* read_indel_error_model_hpp */
