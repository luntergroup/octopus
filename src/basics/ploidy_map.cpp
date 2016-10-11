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
    sample_contigs_[sample].emplace(contig, ploidy);
}

unsigned PloidyMap::organism() const noexcept
{
    return organism_;
}
    
unsigned PloidyMap::operator()(const SampleName& sample, const ContigName& contig) const noexcept
{
    const auto sample_itr = sample_contigs_.find(sample);
    if (sample_itr != std::cend(sample_contigs_)) {
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

} // namespace octopus
