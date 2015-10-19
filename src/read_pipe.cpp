//
//  read_pipe.cpp
//  Octopus
//
//  Created by Daniel Cooke on 19/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "read_pipe.hpp"

#include "read_utils.hpp"

#include <iostream> // DEBUG

namespace Octopus {

ReadPipe::ReadPipe(ReadManager& read_manager, ReadFilterer read_filter,
                   Downsampler downsampler, ReadTransform read_transform)
:
read_manager {read_manager},
read_filter {std::move(read_filter)},
downsampler {std::move(downsampler)},
read_transform {std::move(read_transform)}
{}

ReadMap ReadPipe::fetch_reads(const std::vector<SampleIdType>& samples, const GenomicRegion& region)
{
    auto reads = read_manager.fetch_reads(samples, region);
    
    std::cout << "fetched " << count_reads(reads) << " reads" << std::endl;
    
    auto result = filter_reads(std::move(reads), read_filter).first;
    
    std::cout << "found " << count_reads(result) << " good reads" << std::endl;
    
    // TODO: add downsampler in when it doesn't make so many copies
    //auto result = downsampler(std::move(good_reads));
    
    //std::cout << "downsampled to " << count_reads(result) << " reads" << std::endl;
    
    transform_reads(result, read_transform);
    
    return result;
}

} // namespace Octopus
