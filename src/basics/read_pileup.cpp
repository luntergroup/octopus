// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_pileup.hpp"

#include <iterator>
#include <numeric>

#include "utils/map_utils.hpp"
#include "utils/append.hpp"

namespace octopus {

ReadPileup::ReadPileup(ContigRegion::Position position)
: summaries_ {}
, region_ {position, position + 1}
{
    summaries_.emplace("$", ReadSummaries {});
}

const ContigRegion& ReadPileup::mapped_region() const noexcept
{
    return region_;
}

unsigned ReadPileup::depth() const noexcept
{
    return std::accumulate(std::cbegin(summaries_), std::cend(summaries_), 0u,
                           [] (auto curr, const auto& p) { return curr + p.second.size(); });
}

unsigned ReadPileup::depth(const NucleotideSequence& sequence) const noexcept
{
    return this->summaries(sequence).size();
}

void ReadPileup::add(const AlignedRead& read)
{
    const GenomicRegion region {contig_name(read), region_};
    summaries_[copy_sequence(read, region)].push_back({copy_base_qualities(read, region), read.mapping_quality()});
}

std::vector<ReadPileup::NucleotideSequence> ReadPileup::sequences() const
{
    return extract_keys(summaries_);
}

const ReadPileup::ReadSummaries& ReadPileup::summaries(const NucleotideSequence& sequence) const
{
    const auto itr = summaries_.find(sequence);
    if (itr != std::cend(summaries_)) {
        return itr->second;
    } else {
        return summaries_.at("$");
    }
}

std::vector<ReadPileup::BaseQuality> ReadPileup::base_qualities() const
{
    std::vector<BaseQuality> result {};
    result.reserve(depth());
    for (const auto& p : summaries_) {
        for (const auto& read : p.second) utils::append(read.base_qualities, result);
    }
    return result;
}

std::vector<ReadPileup::BaseQuality> ReadPileup::base_qualities(const NucleotideSequence& sequence) const
{
    std::vector<BaseQuality> result {};
    const auto& summary = this->summaries(sequence);
    result.reserve(summary.size());
    for (const auto& read : summary) utils::append(read.base_qualities, result);
    return result;
}

std::vector<ReadPileup::BaseQuality> ReadPileup::base_qualities_not(const NucleotideSequence& sequence) const
{
    std::vector<BaseQuality> result {};
    result.reserve(depth());
    for (const auto& p : summaries_) {
        if (p.first != sequence && p.first != "$") {
            for (const auto& read : p.second) utils::append(read.base_qualities, result);
        }
    }
    return result;
}

unsigned ReadPileup::sum_base_qualities(const NucleotideSequence& sequence) const
{
    const auto& sequence_summaries = this->summaries(sequence);
    return std::accumulate(std::cbegin(sequence_summaries), std::cend(sequence_summaries), 0u,
                           [] (auto curr, const ReadSummary& summary) {
        return curr + std::accumulate(std::cbegin(summary.base_qualities), std::cend(summary.base_qualities), 0u);
    });
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
