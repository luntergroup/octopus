// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_pileup_hpp
#define read_pileup_hpp

#include <vector>
#include <utility>
#include <functional>

#include "config/common.hpp"
#include "concepts/mappable.hpp"
#include "contig_region.hpp"
#include "genomic_region.hpp"
#include "aligned_read.hpp"

namespace octopus {

class ReadPileup : public Mappable<ReadPileup>
{
public:
    using NucleotideSequence = AlignedRead::NucleotideSequence;
    using BaseQuality        = AlignedRead::BaseQuality;
    using MappingQuality     = AlignedRead::MappingQuality;
    
    struct ReadSummary
    {
        std::vector<BaseQuality> base_qualities;
        MappingQuality mapping_quality;
    };
    using ReadSummaries = std::vector<ReadSummary>;
    
    ReadPileup() = delete;
    
    ReadPileup(ContigRegion::Position position);
    
    ReadPileup(const ReadPileup&)            = default;
    ReadPileup& operator=(const ReadPileup&) = default;
    ReadPileup(ReadPileup&&)                 = default;
    ReadPileup& operator=(ReadPileup&&)      = default;
    
    ~ReadPileup() = default;
    
    const ContigRegion& mapped_region() const noexcept;
    
    void add(const AlignedRead& read);
        
    void summaries(std::function<void(const NucleotideSequence&, const ReadSummaries&)> func) const;

private:
    std::vector<std::pair<NucleotideSequence, ReadSummaries>> summaries_;
    ContigRegion region_;
};

using ReadPileups = std::vector<ReadPileup>;

ReadPileups make_pileups(const ReadContainer& reads, const GenomicRegion& region);

} // namespace octopus

#endif
