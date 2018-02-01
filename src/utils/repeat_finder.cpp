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
find_repeat_regions(std::vector<TandemRepeat>& seeds, const GenomicRegion& region,
                    const InexactRepeatDefinition repeat_definition)
{
    auto itr = std::partition(std::begin(seeds), std::end(seeds),
                              [&] (const auto& repeat) noexcept {
                                  return region_size(repeat) >= repeat_definition.min_exact_repeat_seed_length;
                              });
    itr = std::remove_if(itr, std::end(seeds),
                         [&] (const auto& small_repeat) {
                             return std::none_of(std::begin(seeds), itr, [&] (const auto& seed) {
                                 return std::abs(inner_distance(small_repeat, seed)) <= repeat_definition.max_seed_join_distance;
                             });
                         });
    seeds.erase(itr, std::end(seeds));
    std::sort(std::begin(seeds), std::end(seeds));
    return join(extract_covered_regions(seeds), repeat_definition.max_seed_join_distance);
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
