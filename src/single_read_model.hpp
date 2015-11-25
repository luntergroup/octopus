//
//  single_read_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 25/09/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef single_read_model_hpp
#define single_read_model_hpp

#include <unordered_map>
#include <functional>
#include <cstddef> // size_t

#include "common.hpp"
#include "haplotype.hpp"
#include "aligned_read.hpp"

namespace Octopus
{
    class SingleReadModel
    {
    public:
        SingleReadModel()  = default;
        ~SingleReadModel() = default;
        
        SingleReadModel(const SingleReadModel&)            = default;
        SingleReadModel& operator=(const SingleReadModel&) = default;
        SingleReadModel(SingleReadModel&&)                 = default;
        SingleReadModel& operator=(SingleReadModel&&)      = default;
        
        double log_probability(const AlignedRead& read, const Haplotype& haplotype); // ln p(read | haplotype)
    };
    
} // namespace Octopus

#endif /* single_read_model_hpp */
