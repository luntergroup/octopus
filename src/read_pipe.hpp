//
//  read_pipe.hpp
//  Octopus
//
//  Created by Daniel Cooke on 19/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef read_pipe_hpp
#define read_pipe_hpp

#include <vector>
#include <unordered_map>
#include <cstddef>
#include <functional>

#include <boost/optional.hpp>

#include "common.hpp"
#include "genomic_region.hpp"
#include "read_manager.hpp"
#include "read_filterer.hpp"
#include "read_transformer.hpp"
#include "downsampler.hpp"
#include "logging.hpp"

namespace octopus
{
/*
 ReadPipe provides a wrapper for the basic manipulation of AlignedRead's. It is responsable for
 actually fetching the reads from file, and applying any filters, transforms etc. The result is
 then a set of reads that can be used for calling.
 
 Rather than just fetching all reads for all samples in one go and processing, we can be a bit more
 clever and fetch reads in batches (e.g. each sample, or regions in samples). This could potentially
 decrease average memory consumption (and also increase runtime performance) by minimising the
 number of 'bad' reads in memory. If we are really short on memory we could even compress filtered
 read batches while we process other batches.
 */
class ReadPipe
{
public:
    using ReadTransformer = preprocess::transform::ReadTransformer;
    using ReadFilterer    = preprocess::filter::ReadFiltererTp<ReadManager::ReadContainer>;
    using Downsampler     = preprocess::Downsampler;
    
    //struct Report;
    
    ReadPipe() = delete;
    
    ReadPipe(const ReadManager& manager, std::vector<SampleName> samples);
    
    ReadPipe(const ReadManager& manager, ReadTransformer transformer,
             ReadFilterer filterer, boost::optional<Downsampler> downsampler,
             std::vector<SampleName> samples);
    
    ReadPipe(const ReadPipe&)            = delete;
    ReadPipe& operator=(const ReadPipe&) = delete;
    ReadPipe(ReadPipe&&)                 = default;
    ReadPipe& operator=(ReadPipe&&)      = default;
    
    ~ReadPipe() = default;
    
    void set_read_manager(const ReadManager& manager) noexcept;
    
    unsigned num_samples() const noexcept;
    const std::vector<SampleName>& samples() const noexcept;
    
    ReadMap fetch_reads(const GenomicRegion& region) const;
    ReadMap fetch_reads(const std::vector<GenomicRegion>& regions) const;
    
    //Report get_report() const;
    
private:
    std::reference_wrapper<const ReadManager> manager_;
    ReadTransformer transformer_;
    ReadFilterer filterer_;
    boost::optional<Downsampler> downsampler_;
    
    std::vector<SampleName> samples_;
    
    mutable boost::optional<logging::DebugLogger> debug_log_;
};
} // namespace octopus

#endif /* read_pipe_hpp */
