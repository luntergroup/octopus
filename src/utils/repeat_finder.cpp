// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "repeat_finder.hpp"

namespace octopus {

std::vector<TandemRepeat>
find_exact_tandem_repeats(const ReferenceGenome& reference, const GenomicRegion& region, unsigned max_period)
{
    auto sequence = reference.fetch_sequence(region);
    return find_exact_tandem_repeats(sequence, region, 1, max_period);
}

std::vector<GenomicRegion>
find_repeat_regions(const std::vector<TandemRepeat>& repeats, const GenomicRegion& region,
                    const InexactRepeatDefinition repeat_definition)
{
    assert(std::is_sorted(std::cbegin(repeats), std::cend(repeats)));
    std::vector<TandemRepeat> seeds {};
    seeds.reserve(repeats.size());
    std::copy_if(std::cbegin(repeats), std::cend(repeats), std::back_inserter(seeds),
                 [&] (const auto& repeat) { return region_size(repeat) >= repeat_definition.min_exact_repeat_seed_length; });
    auto repeat_begin_itr = std::cbegin(repeats);
    std::vector<GenomicRegion> hits {};
    hits.reserve(repeats.size());
    for (const auto& seed : seeds) {
        const auto expanded_seed_region = expand(seed.region, repeat_definition.max_seed_join_distance);
        const auto overlapped_repeats = overlap_range(repeat_begin_itr, std::cend(repeats), expanded_seed_region);
        for (const auto& repeat : overlapped_repeats) {
            hits.push_back(repeat.region);
        }
        repeat_begin_itr = overlapped_repeats.begin().base();
    }
    auto result = join(extract_covered_regions(hits), repeat_definition.max_seed_join_distance);
    result.erase(std::remove_if(std::begin(result), std::end(result),
                                [&] (const auto& region) {
                                    return size(region) < repeat_definition.min_joined_repeat_length;
                                }), std::end(result));
    return result;
}

std::vector<GenomicRegion>
find_repeat_regions(const ReferenceGenome& reference, const GenomicRegion& region,
                    const InexactRepeatDefinition repeat_definition)
{
    auto sequence = reference.fetch_sequence(region);
    auto seeds = find_exact_tandem_repeats(sequence, region, 1, repeat_definition.max_exact_repeat_seed_period);
    return find_repeat_regions(seeds, region, repeat_definition);
}

} // namespace octopus
