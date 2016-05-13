//
//  read_utils.cpp
//  Octopus
//
//  Created by Daniel Cooke on 09/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "read_utils.hpp"

namespace Octopus
{

// TODO
AlignedRead find_next_segment(const AlignedRead& read,
                              const MappableMap<GenomicRegion::ContigNameType, AlignedRead>& reads)
{
    if (!read.is_chimeric()) {
        throw std::runtime_error {"cannot find next segment as read is not chimeric"};
    }
    
    auto segment_region = read.next_segment().inferred_region();
    
    return read;
}

// TODO
MappableFlatMultiSet<AlignedRead> find_chimeras(const AlignedRead& read, const MappableFlatMultiSet<AlignedRead>& reads)
{
    return MappableFlatMultiSet<AlignedRead> {};
}

} // namespace Octopus
