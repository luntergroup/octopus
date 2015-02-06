//
//  genome_region.cpp
//  Octopus
//
//  Created by Daniel Cooke on 06/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "genome_region.h"

GenomeRegion::GenomeRegion(std::string contig_name, size_t begin, size_t end) noexcept
: contig_name {contig_name}, begin {begin}, end {end}
{}
