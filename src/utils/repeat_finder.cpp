// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "repeat_finder.hpp"

namespace octopus {

std::vector<TandemRepeat>
find_exact_tandem_repeats(const ReferenceGenome& reference, const GenomicRegion& region, unsigned max_period)
{
    auto sequence = reference.fetch_sequence(region);
    return find_exact_tandem_repeats(sequence, region, 1, max_period);
}

bool is_good_seed(const TandemRepeat& repeat, const InexactRepeatDefinition& repeat_def) noexcept
{
    const auto repeat_length = region_size(repeat);
    return repeat_length >= repeat_def.min_exact_repeat_seed_length
           && repeat_length / repeat.period() >= repeat_def.min_exact_repeat_seed_periods;
}

std::vector<GenomicRegion>
find_repeat_regions(const std::vector<TandemRepeat>& repeats, const GenomicRegion& region,
                    const InexactRepeatDefinition repeat_def)
{
    assert(std::is_sorted(std::cbegin(repeats), std::cend(repeats)));
    std::vector<TandemRepeat> seeds {};
    seeds.reserve(repeats.size());
    std::copy_if(std::cbegin(repeats), std::cend(repeats), std::back_inserter(seeds),
                 [&] (const auto& repeat) { return is_good_seed(repeat, repeat_def); });
    if (seeds.size() < repeat_def.min_exact_seeds) {
        return {};
    }
    auto repeat_begin_itr = std::cbegin(repeats);
    std::vector<GenomicRegion> hits {};
    hits.reserve(repeats.size());
    for (const auto& seed : seeds) {
        const auto expanded_seed_region = expand(mapped_region(seed), repeat_def.max_seed_join_distance);
        const auto overlapped_repeats = overlap_range(repeat_begin_itr, std::cend(repeats), expanded_seed_region);
        for (const auto& repeat : overlapped_repeats) {
            hits.push_back(mapped_region(repeat));
        }
        repeat_begin_itr = overlapped_repeats.begin().base();
    }
    auto result = join(extract_covered_regions(hits), repeat_def.max_seed_join_distance);
    result.erase(std::remove_if(std::begin(result), std::end(result),
                                [&] (const auto& region) {
                                    return size(region) < repeat_def.min_joined_repeat_length;
                                }), std::end(result));
    return result;
}

std::vector<GenomicRegion>
find_repeat_regions(const ReferenceGenome& reference, const GenomicRegion& region,
                    const InexactRepeatDefinition repeat_def)
{
    auto sequence = reference.fetch_sequence(region);
    auto seeds = find_exact_tandem_repeats(sequence, region, 1, repeat_def.max_exact_repeat_seed_period);
    return find_repeat_regions(seeds, region, repeat_def);
}

} // namespace octopus
