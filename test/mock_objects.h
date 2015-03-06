//
//  mock_objects.h
//  Octopus
//
//  Created by Daniel Cooke on 22/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_mock_objects_h
#define Octopus_mock_objects_h

#include "genomic_region.h"
#include "aligned_read.h"
#include "cigar_string.h"

inline GenomicRegion get_mock_region()
{
    return GenomicRegion {"1", 0, 1};
}

inline AlignedRead get_mock_aligned_read(std::string the_sequence)
{
    return AlignedRead {
        get_mock_region(),
        std::move(the_sequence),
        AlignedRead::Qualities {},
        CigarString {},
        0,
        "1",
        0,
        0
    };
}

#endif
