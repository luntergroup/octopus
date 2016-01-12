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
    
    auto segment_region = read.get_next_segment().get_inferred_region();
    
    return read;
}

// TODO
MappableSet<AlignedRead> find_chimeras(const AlignedRead& read, const MappableSet<AlignedRead>& reads)
{
    return MappableSet<AlignedRead> {};
}

} // namespace Octopus
