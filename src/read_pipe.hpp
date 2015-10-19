//
//  read_pipe.hpp
//  Octopus
//
//  Created by Daniel Cooke on 19/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef read_pipe_hpp
#define read_pipe_hpp

#include "common.hpp"
#include "read_manager.hpp"
#include "read_filter.hpp"
#include "read_transform.hpp"
#include "downsampler.hpp"

namespace Octopus {

class ReadPipe
{
public:
    ReadPipe() = delete;
    ReadPipe(ReadManager& read_manager, ReadFilterer read_filter,
             Downsampler downsampler, ReadTransform read_transform);
    
    ReadMap fetch_reads(const std::vector<SampleIdType>& samples, const GenomicRegion& region);
    
private:
    ReadManager& read_manager;
    ReadFilterer read_filter;
    Downsampler downsampler;
    ReadTransform read_transform;
};
    
} // namespace Octopus

#endif /* read_pipe_hpp */
