// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_pipe_hpp
#define read_pipe_hpp

#include <vector>
#include <unordered_map>
#include <cstddef>
#include <functional>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "basics/genomic_region.hpp"
#include "io/read/read_manager.hpp"
#include "logging/logging.hpp"
#include "filtering/read_filterer.hpp"
#include "transformers/read_transformer.hpp"
#include "downsampling/downsampler.hpp"

namespace octopus {
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
    using ReadTransformer = readpipe::ReadTransformer;
    using ReadFilterer    = readpipe::ReadFiltererTp<ReadManager::ReadContainer>;
    using Downsampler     = readpipe::Downsampler;
    
    struct Report
    {
        readpipe::DownsamplerReportMap downsample_report;
    };
    
    ReadPipe() = delete;
    
    ReadPipe(const ReadManager& source, std::vector<SampleName> samples);
    
    ReadPipe(const ReadManager& source, ReadTransformer transformer,
             ReadFilterer filterer, boost::optional<Downsampler> downsampler,
             std::vector<SampleName> samples);
    
    ReadPipe(const ReadManager& manager, ReadTransformer prefilter_transformer,
             ReadFilterer filterer, ReadTransformer postfilter_transformer,
             boost::optional<Downsampler> downsampler,
             std::vector<SampleName> samples);
    
    ReadPipe(const ReadPipe&)            = delete;
    ReadPipe& operator=(const ReadPipe&) = delete;
    ReadPipe(ReadPipe&&)                 = default;
    ReadPipe& operator=(ReadPipe&&)      = default;
    
    ~ReadPipe() = default;
    
    const ReadManager& read_manager() const noexcept;
    void set_read_manager(const ReadManager& source) noexcept;
    
    unsigned num_samples() const noexcept;
    const std::vector<SampleName>& samples() const noexcept;
    
    ReadMap fetch_reads(const GenomicRegion& region, boost::optional<Report&> report = boost::none) const;
    ReadMap fetch_reads(const std::vector<GenomicRegion>& regions, boost::optional<Report&> report = boost::none) const;
    
    //Report get_report() const;
    
private:
    std::reference_wrapper<const ReadManager> source_;
    ReadTransformer prefilter_transformer_;
    ReadFilterer filterer_;
    boost::optional<ReadTransformer> postfilter_transformer_;
    boost::optional<Downsampler> downsampler_;
    std::vector<SampleName> samples_;
    mutable boost::optional<logging::DebugLogger> debug_log_;
};

} // namespace octopus

#endif
