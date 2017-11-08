// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef repeat_finer_hpp
#define repeat_finer_hpp

#include <vector>
#include <iterator>
#include <algorithm>
#include <cmath>

#include "io/reference/reference_genome.hpp"
#include "basics/genomic_region.hpp"
#include "concepts/mappable.hpp"
#include "utils/mappable_algorithms.hpp"
#include "tandem/tandem.hpp"

namespace octopus {

struct TandemRepeat : public Mappable<TandemRepeat>
{
    using SizeType = GenomicRegion::Size;
    TandemRepeat() = delete;
    template <typename T>
    TandemRepeat(T region, GenomicRegion::Size period)
    : region {std::forward<T>(region)}, period {period} {}
    
    GenomicRegion region;
    GenomicRegion::Size period;
    
    const GenomicRegion& mapped_region() const noexcept { return region; }
};

struct InexactRepeatDefinition
{
    unsigned max_exact_repeat_seed_period = 100;
    unsigned min_exact_repeat_seed_length = 100;
    unsigned max_seed_join_distance       = 50;
};

template <typename SequenceType>
std::vector<TandemRepeat>
find_exact_tandem_repeats(SequenceType sequence, const GenomicRegion& region,
                          GenomicRegion::Size min_period = 2, GenomicRegion::Size max_period = 10000)
{
    if (sequence.back() != 'N') {
        sequence.reserve(sequence.size() + 1);
        sequence.push_back('N');
    }
    auto n_shift_map = tandem::collapse(sequence, 'N');
    auto maximal_repetitions = tandem::extract_exact_tandem_repeats(sequence , min_period, max_period);
    tandem::rebase(maximal_repetitions, n_shift_map);
    n_shift_map.clear();
    std::vector<TandemRepeat> result {};
    result.reserve(maximal_repetitions.size());
    auto offset = region.begin();
    for (const auto& run : maximal_repetitions) {
        result.emplace_back(GenomicRegion {region.contig_name(),
                                           static_cast<GenomicRegion::Size>(run.pos + offset),
                                           static_cast<GenomicRegion::Size>(run.pos + run.length + offset)
        }, run.period);
    }
    return result;
}

template <typename SequenceType>
std::vector<GenomicRegion>
find_repeat_regions(const SequenceType& sequence, const GenomicRegion& region,
                    const InexactRepeatDefinition repeat_definition)
{
    auto repeats = find_exact_tandem_repeats(sequence, region, 2, repeat_definition.max_exact_repeat_seed_period);
    auto itr = std::partition(std::begin(repeats), std::end(repeats),
                              [&] (const auto& repeat) noexcept {
                                  return region_size(repeat) >= repeat_definition.min_exact_repeat_seed_length;
                              });
    itr = std::remove_if(itr, std::end(repeats),
                         [&] (const auto& small_repeat) {
                             return std::none_of(std::begin(repeats), itr, [&] (const auto& seed) {
                                 return std::abs(inner_distance(small_repeat, seed)) <= repeat_definition.max_seed_join_distance;
                             });
                         });
    repeats.erase(itr, std::end(repeats));
    std::sort(std::begin(repeats), std::end(repeats));
    return join(extract_covered_regions(repeats), repeat_definition.max_seed_join_distance);
}

std::vector<GenomicRegion>
find_repeat_regions(const ReferenceGenome& reference, const GenomicRegion& region,
                    InexactRepeatDefinition repeat_definition = InexactRepeatDefinition {});

} // namespace octopus

#endif
