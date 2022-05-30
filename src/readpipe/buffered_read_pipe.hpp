// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef buffered_read_pipe_hpp
#define buffered_read_pipe_hpp

#include <functional>
#include <cstddef>

#include <boost/optional.hpp>

#include "read_pipe.hpp"
#include "basics/genomic_region.hpp"
#include "containers/mappable_map.hpp"
#include "logging/logging.hpp"

namespace octopus {

class BufferedReadPipe
{
public:
    struct Config
    {
        std::size_t max_buffer_size;
        GenomicRegion::Size fetch_expansion = 0;
        boost::optional<GenomicRegion::Size> max_fetch_size = boost::none;
        boost::optional<GenomicRegion::Size> max_hint_gap = boost::none;
        bool allow_unchecked_fetches = true;
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
    ReadMap fetch_reads(const std::vector<GenomicRegion>& regions) const;
    
    // Hints are regions that will likely be requested in the future. This does not affect the observable behaviour of
    // the object, but may allow improved performance through optimised read buffering. If the hints given are
    // inaccurate it will likely result in worse performance.
    void hint(std::vector<GenomicRegion> hints) const;
    
    bool is_cached(const GenomicRegion& region) const noexcept;
    
private:
    using RegionMap = MappableSetMap<GenomicRegion::ContigName, GenomicRegion>;
    
    std::reference_wrapper<const ReadPipe> source_;
    Config config_;
    mutable ReadMap buffer_;
    mutable boost::optional<GenomicRegion> buffered_region_;
    mutable RegionMap hints_;
    mutable bool default_unchecked_fetch_overflowed_ = false;
    mutable bool adjusted_unchecked_fetch_overflowed_ = false;
    mutable boost::optional<GenomicRegion::Size> min_checked_fetch_size_ = boost::none;
    mutable boost::optional<logging::DebugLogger> debug_log_;
    
    void setup_buffer(const GenomicRegion& request) const;
    GenomicRegion get_max_fetch_region(const GenomicRegion& request) const;
    GenomicRegion get_default_max_fetch_region(const GenomicRegion& request) const;
    bool can_make_unchecked_fetch() const noexcept;
};

} // namespace octopus

#endif
