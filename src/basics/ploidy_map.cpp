// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "ploidy_map.hpp"

#include <iterator>
#include <limits>

namespace octopus {

PloidyMap::PloidyMap(unsigned organism)
: organism_ {organism}
{}

void PloidyMap::set(const ContigName& contig, unsigned ploidy)
{
    contigs_.emplace(contig, ploidy);
}

void PloidyMap::set(const SampleName& sample, const ContigName& contig, unsigned ploidy)
{
    allosomes_[contig].emplace(sample, ploidy);
}

bool PloidyMap::is_autosome(const ContigName& contig) const noexcept
{
    return allosomes_.count(contig) == 0;
}

unsigned PloidyMap::of(const SampleName& sample, const ContigName& contig) const noexcept
{
    const auto sample_itr = allosomes_.find(contig);
    if (sample_itr != std::cend(allosomes_)) {
        const auto contig_itr = sample_itr->second.find(sample);
        if (contig_itr != std::cend(sample_itr->second)) {
            return contig_itr->second;
        }
    }
    const auto contig_itr = contigs_.find(contig);
    if (contig_itr != std::cend(contigs_)) {
        return contig_itr->second;
    }
    return organism_;
}

std::vector<unsigned> get_ploidies(const std::vector<SampleName>& samples, const ContigName& contig, const PloidyMap& ploidies)
{
    std::vector<unsigned> result(samples.size());
    if (ploidies.is_autosome(contig)) {
        result.assign(samples.size(), ploidies.of(samples.front(), contig));
    } else {
        std::transform(std::cbegin(samples), std::cend(samples), std::begin(result),
                       [&] (const auto& sample) { return ploidies.of(sample, contig); });
    }
    return result;
}

unsigned get_min_ploidy(const std::vector<SampleName>& samples, const std::vector<ContigName>& contigs, const PloidyMap& ploidies)
{
    auto result = std::numeric_limits<unsigned>::max();
    for (const auto& sample : samples) {
        for (const auto& contig : contigs) {
            result = std::min(result, ploidies.of(sample, contig));
        }
    }
    return result;
}

unsigned get_max_ploidy(const std::vector<SampleName>& samples, const std::vector<ContigName>& contigs, const PloidyMap& ploidies)
{
    unsigned result {0};
    for (const auto& sample : samples) {
        for (const auto& contig : contigs) {
            result = std::max(result, ploidies.of(sample, contig));
        }
    }
    return result;
}

} // namespace octopus
