// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef buffered_read_pipe_hpp
#define buffered_read_pipe_hpp

#include <functional>
#include <cstddef>

#include <boost/optional.hpp>

#include "read_pipe.hpp"
#include "basics/genomic_region.hpp"

namespace octopus {

class BufferedReadPipe
{
public:
    BufferedReadPipe() = delete;
    
    BufferedReadPipe(ReadPipe& source, std::size_t max_buffer_size);
    
    BufferedReadPipe(const BufferedReadPipe&)            = delete;
    BufferedReadPipe& operator=(const BufferedReadPipe&) = delete;
    BufferedReadPipe(BufferedReadPipe&&)                 = default;
    BufferedReadPipe& operator=(BufferedReadPipe&&)      = default;
    
    ~BufferedReadPipe() = default;
    
    ReadPipe& source() noexcept;
    
    void clear() noexcept;
    
    ReadMap fetch_reads(const GenomicRegion& region) const;
    
private:
    std::reference_wrapper<ReadPipe> source_;
    std::size_t max_buffer_size_;
    mutable ReadMap buffer_;
    mutable boost::optional<GenomicRegion> buffered_region_;
    
    bool requires_new_fetch(const GenomicRegion& region) const noexcept ;
    void setup_buffer(const GenomicRegion& region) const;
};

} // namespace octopus

#endif
