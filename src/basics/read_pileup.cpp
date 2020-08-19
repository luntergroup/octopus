// Copyright (c) 2015-2020 Daniel Cooke
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
    auto sequence = copy_sequence(read, region);
    auto itr = std::find_if(std::next(std::begin(summaries_)), std::end(summaries_), [&] (const auto& p) { return p.first == sequence; });
    if (itr == std::cend(summaries_)) {
        summaries_.emplace_back(std::move(sequence), ReadSummaries {});
        itr = std::prev(std::end(summaries_));
    }
    itr->second.push_back({copy_base_qualities(read, region), read.mapping_quality()});
}

std::vector<ReadPileup::NucleotideSequence> ReadPileup::sequences() const
{
    std::vector<NucleotideSequence> result {};
    result.reserve(summaries_.size() - 1);
    std::transform(std::next(std::cbegin(summaries_)), std::cend(summaries_), std::back_inserter(result),
                   [] (const auto& p) { return p.first; });
    return result;
}

const ReadPileup::ReadSummaries& ReadPileup::summaries(const NucleotideSequence& sequence) const
{
    const auto itr = std::find_if(std::next(std::cbegin(summaries_)), std::cend(summaries_), [&] (const auto& p) { return p.first == sequence; });
    if (itr != std::cend(summaries_)) {
        return itr->second;
    } else {
        return summaries_.front().second;
    }
}

std::vector<ReadPileup::BaseQuality> ReadPileup::base_qualities() const
{
    std::vector<BaseQuality> result {};
    if (summaries_.size() == 1) return result;
    result.reserve(depth());
    for (const auto& p : summaries_) {
        for (const auto& read : p.second) utils::append(read.base_qualities, result);
    }
    return result;
}

std::vector<ReadPileup::BaseQuality> ReadPileup::base_qualities(const NucleotideSequence& sequence) const
{
    std::vector<BaseQuality> result {};
    if (summaries_.size() == 1) return result;
    const auto& summary = this->summaries(sequence);
    result.reserve(summary.size());
    for (const auto& read : summary) utils::append(read.base_qualities, result);
    return result;
}

std::vector<ReadPileup::BaseQuality> ReadPileup::base_qualities_not(const NucleotideSequence& sequence) const
{
    std::vector<BaseQuality> result {};
    if (summaries_.size() == 1) return result;
    result.reserve(depth());
    std::for_each(std::next(std::cbegin(summaries_)), std::cend(summaries_), [&] (const auto& p) {
        if (p.first != sequence) {
            for (const auto& read : p.second) utils::append(read.base_qualities, result);
        }
    });
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
