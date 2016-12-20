// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "ploidy_map.hpp"

#include <iterator>

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
    allosomes_[sample].emplace(contig, ploidy);
}

bool PloidyMap::is_autosome(const ContigName& contig) const noexcept
{
    return allosomes_.count(contig) == 0;
}
    
unsigned PloidyMap::of(const SampleName& sample, const ContigName& contig) const noexcept
{
    const auto sample_itr = allosomes_.find(sample);
    if (sample_itr != std::cend(allosomes_)) {
        const auto contig_itr = sample_itr->second.find(contig);
        if (contig_itr != std::cend(sample_itr->second)) {
            return contig_itr->second;
        }
    } else {
        const auto contig_itr = contigs_.find(contig);
        if (contig_itr != std::cend(contigs_)) {
            return contig_itr->second;
        }
    }
    return organism_;
}

std::vector<unsigned> get_ploidies(const std::vector<SampleName>& samples, const ContigName& contig, const PloidyMap& ploidies)
{
    std::vector<unsigned> result {};
    if (ploidies.is_autosome(contig)) {
        const auto ploidy = ploidies.of(samples.front(), contig);
        result.assign(samples.size(), ploidy);
    } else {
        result.reserve(samples.size());
        std::transform(std::cbegin(samples), std::cend(samples), std::back_inserter(result),
                       [&] (const auto& sample) { return ploidies.of(sample, contig); });
    }
    return result;
}

} // namespace octopus
