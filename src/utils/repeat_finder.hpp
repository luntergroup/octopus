// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef repeat_finder_hpp
#define repeat_finder_hpp

#include <vector>
#include <iterator>
#include <algorithm>
#include <cmath>

#include "io/reference/reference_genome.hpp"
#include "basics/genomic_region.hpp"
#include "basics/tandem_repeat.hpp"
#include "concepts/mappable.hpp"
#include "concepts/comparable.hpp"
#include "utils/mappable_algorithms.hpp"
#include "tandem/tandem.hpp"

namespace octopus {

struct InexactRepeatDefinition
{
    unsigned max_exact_repeat_seed_period  = 6;
    unsigned min_exact_repeat_seed_length  = 3;
    unsigned min_exact_repeat_seed_periods = 3;
    unsigned min_exact_seeds               = 1;
    unsigned max_seed_join_distance        = 2;
    unsigned min_joined_repeat_length      = 10;
};

template <typename SequenceType>
std::vector<TandemRepeat>
find_exact_tandem_repeats(SequenceType& sequence, const GenomicRegion& region,
                          GenomicRegion::Size min_period, GenomicRegion::Size max_period)
{
    if (sequence.back() != '$') {
        sequence.reserve(sequence.size() + 1);
        sequence.push_back('$');
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
        }, SequenceType {std::next(std::cbegin(sequence), run.pos), std::next(std::cbegin(sequence), run.pos + run.period)});
    }
    return result;
}

template <typename SequenceType>
std::vector<TandemRepeat>
find_exact_tandem_repeats(const SequenceType& sequence, const GenomicRegion& region,
                          GenomicRegion::Size min_period, GenomicRegion::Size max_period)
{
    auto tmp = sequence;
    return find_exact_tandem_repeats(tmp, region, min_period, max_period);
}

std::vector<TandemRepeat>
find_exact_tandem_repeats(const ReferenceGenome& reference, const GenomicRegion& region, unsigned max_period);

std::vector<GenomicRegion>
find_repeat_regions(const std::vector<TandemRepeat>& repeats, const GenomicRegion& region,
                    const InexactRepeatDefinition repeat_def);

std::vector<GenomicRegion>
find_repeat_regions(const ReferenceGenome& reference, const GenomicRegion& region,
                    InexactRepeatDefinition repeat_def = InexactRepeatDefinition {});

} // namespace octopus

#endif
