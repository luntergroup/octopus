// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_pileup.hpp"

namespace octopus {

auto make_pileup(const ReadContainer& reads, const GenomicRegion& region)
{
    ReadPileup result {region};
    for (const auto& read : overlap_range(reads, region)) {
        auto subsequence = copy_sequence(read, region);
        auto base_qualities = copy_base_qualities(read, region);
        result.read_sequences[subsequence].push_back({std::move(base_qualities), read.mapping_quality()});
    }
    return result;
}

} // namespace octopus
