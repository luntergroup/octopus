// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_pileup.hpp"

#include <iterator>
#include <numeric>

namespace octopus {

ReadPileup::ReadPileup(ContigRegion::Position position) : region_ {position, position + 1} {}

const ContigRegion& ReadPileup::mapped_region() const noexcept
{
    return region_;
}

unsigned ReadPileup::depth() const noexcept
{
    return std::accumulate(std::cbegin(summaries_), std::cend(summaries_), 0u,
                           [] (auto curr, const auto& p) { return curr + p.second.size(); });
}

void ReadPileup::add(const AlignedRead& read)
{
    const GenomicRegion region {contig_name(read), region_};
    summaries_[copy_sequence(read, region)].push_back({copy_base_qualities(read, region), read.mapping_quality()});
}

namespace {

auto overlap_range(std::vector<ReadPileup>& pileups, const AlignedRead& read)
{
    return overlap_range(std::begin(pileups), std::end(pileups), contig_region(read), BidirectionallySortedTag {});
}

} // namespace

ReadPileups make_pileups(const ReadContainer& reads, const GenomicRegion& region)
{
    ReadPileups result {};
    result.reserve(size(region));
    for (auto position = region.begin(); position < region.end(); ++position) {
        result.emplace_back(position);
    }
    for (const AlignedRead& read : overlap_range(reads, region)) {
        for (ReadPileup& pileup : overlap_range(result, read)) {
            pileup.add(read);
        }
    }
    return result;
}

} // namespace octopus
