// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef buffered_read_pipe_hpp
#define buffered_read_pipe_hpp

#include <functional>
#include <cstddef>

#include <boost/optional.hpp>

#include "read_pipe.hpp"
#include "basics/genomic_region.hpp"
#include "containers/mappable_map.hpp"

namespace octopus {

class BufferedReadPipe
{
public:
    struct Config
    {
        std::size_t max_buffer_size;
        boost::optional<GenomicRegion::Size> max_fetch_size = boost::none;
        boost::optional<GenomicRegion::Size> max_hint_gap = boost::none;
    };
    
    BufferedReadPipe() = delete;
    
    BufferedReadPipe(const ReadPipe& source, Config config);
    BufferedReadPipe(const ReadPipe& source, Config config, std::vector<GenomicRegion> hints);
    
    BufferedReadPipe(const BufferedReadPipe&)            = delete;
    BufferedReadPipe& operator=(const BufferedReadPipe&) = delete;
    BufferedReadPipe(BufferedReadPipe&&)                 = default;
    BufferedReadPipe& operator=(BufferedReadPipe&&)      = default;
    
    ~BufferedReadPipe() = default;
    
    const ReadPipe& source() const noexcept;
    
    void clear() noexcept;
    
    ReadMap fetch_reads(const GenomicRegion& region) const;
    
    // Hints are regions that will likely be requested in the future. This does not affect the observable behaviour of
    // the object, but may allow improved performance through optimised read buffering. If the hints given are
    // inaccurate it will likely result in worse performance.
    void hint(std::vector<GenomicRegion> hints) const;
    
private:
    using RegionMap = MappableSetMap<GenomicRegion::ContigName, GenomicRegion>;
    
    std::reference_wrapper<const ReadPipe> source_;
    Config config_;
    mutable ReadMap buffer_;
    mutable boost::optional<GenomicRegion> buffered_region_;
    mutable RegionMap hints_;
    
    bool requires_new_fetch(const GenomicRegion& request) const noexcept ;
    void setup_buffer(const GenomicRegion& request) const;
    GenomicRegion get_max_fetch_region(const GenomicRegion& request) const;
    GenomicRegion get_default_max_fetch_region(const GenomicRegion& request) const;
};

} // namespace octopus

#endif
