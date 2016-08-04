//
//  trio_caller.hpp
//  octopus
//
//  Created by Daniel Cooke on 03/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef trio_caller_hpp
#define trio_caller_hpp

#include <vector>
#include <string>

#include "caller.hpp"

namespace octopus {
    
    class GenomicRegion;
    class Variant;
    class VcfRecord;
    class ReadPipe;
    
    class TrioCaller// : public Caller
    {
    public:
        TrioCaller() = delete;
        
        ~TrioCaller() = default;
        
        TrioCaller(const TrioCaller&)            = delete;
        TrioCaller& operator=(const TrioCaller&) = delete;
        TrioCaller(TrioCaller&&)                 = delete;
        TrioCaller& operator=(TrioCaller&&)      = delete;
        
    private:
        const unsigned ploidy_;
        const SampleName mother_, father_;
        const double min_variant_posterior_ = 0.95;
    };
    
} // namespace octopus

#endif /* trio_caller_hpp */
