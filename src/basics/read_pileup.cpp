// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_pileup.hpp"

#include <iterator>
#include <algorithm>
#include <numeric>

#include "utils/append.hpp"
#include "utils/mappable_algorithms.hpp"

namespace octopus {

ReadPileup::ReadPileup(ContigRegion::Position position)
: summaries_ {}
, region_ {position, position + 1}
{
    summaries_.reserve(3);
    summaries_.emplace_back("$", ReadSummaries {});
}

const ContigRegion& ReadPileup::mapped_region() const noexcept
{
    return region_;
}

void ReadPileup::add(const AlignedRead& read)
{
    const GenomicRegion region {contig_name(read), region_};
    auto sequence = copy_sequence(read, region);
    auto itr = std::find_if(std::next(std::begin(summaries_)), std::end(summaries_), [&] (const auto& p) { return p.first == sequence; });
    if (itr == std::cend(summaries_)) {
        summaries_.emplace_back(std::move(sequence), ReadSummaries {});
        itr = std::prev(std::end(summaries_));
    }
    itr->second.push_back({copy_base_qualities(read, region), read.mapping_quality()});
}

void ReadPileup::summaries(std::function<void(const NucleotideSequence&, const ReadSummaries&)> func) const
{
    for (const auto& p : summaries_) {
        func(p.first, p.second);
    }
}

namespace {

auto pileup_overlap_range(std::vector<ReadPileup>& pileups, const AlignedRead& read)
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
        for (ReadPileup& pileup : pileup_overlap_range(result, read)) {
            pileup.add(read);
        }
    }
    return result;
}

} // namespace octopus
