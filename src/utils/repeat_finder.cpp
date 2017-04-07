// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "repeat_finder.hpp"

namespace octopus {

std::vector<GenomicRegion>
find_repeat_regions(const ReferenceGenome& reference, const GenomicRegion& region,
                    const InexactRepeatDefinition repeat_definition)
{
    const auto sequence = reference.fetch_sequence(region);
    return find_repeat_regions(sequence, region, repeat_definition);
}

} // namespace octopus
