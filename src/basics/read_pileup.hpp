// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_pileup_hpp
#define read_pileup_hpp

#include <vector>
#include <map>
#include <utility>

#include "config/common.hpp"
#include "concepts/mappable.hpp"
#include "genomic_region.hpp"
#include "aligned_read.hpp"

namespace octopus {

struct ReadPileup : public Mappable<ReadPileup>
{
    using NucleotideSequence = AlignedRead::NucleotideSequence;
    using BaseQuality        = AlignedRead::BaseQuality;
    using MappingQuality     = AlignedRead::MappingQuality;
    struct ReadSummary
    {
        std::vector<BaseQuality> base_qualities;
        MappingQuality mapping_quality;
    };
    std::map<NucleotideSequence, std::vector<ReadSummary>> read_sequences;
    GenomicRegion region;
    ReadPileup(GenomicRegion region) : region { std::move(region) } {}
    const GenomicRegion& mapped_region() const noexcept { return region; }
};

auto make_pileup(const ReadContainer& reads, const GenomicRegion& region);

} // namespace octopus

#endif
