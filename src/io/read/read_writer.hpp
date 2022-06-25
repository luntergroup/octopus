// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_writer_hpp
#define read_writer_hpp

#include <memory>
#include <mutex>

#include <boost/filesystem/path.hpp>

#include "concepts/equitable.hpp"
#include "htslib_sam_facade.hpp"

namespace octopus {

class AlignedRead;

namespace io {

class ReadWriter
{
public:
    using Path = boost::filesystem::path;
    
    ReadWriter() = delete;
    
    ReadWriter(Path bam_out, Path bam_template);
    
    ReadWriter(const ReadWriter&)            = delete;
    ReadWriter& operator=(const ReadWriter&) = delete;
    ReadWriter(ReadWriter&&);
    ReadWriter& operator=(ReadWriter&&)      = delete;
    
    ~ReadWriter() = default;
    
    friend void swap(ReadWriter& lhs, ReadWriter& rhs) noexcept;
    
    const Path& path() const noexcept;
    
    void write(const AlignedRead& read);
    
private:
    Path path_;
    std::unique_ptr<HtslibSamFacade> impl_;
    mutable std::mutex mutex_;
};

ReadWriter& operator<<(ReadWriter& dst, const AlignedRead& read);

template <typename Container>
void write(const Container& reads, ReadWriter& dst)
{
    for (const auto& read : reads) {
        dst << read;
    }
}

template <typename Container>
ReadWriter& operator<<(ReadWriter& dst, const Container& reads)
{
    write(reads, dst);
    return dst;
}

} // namespace io
} // namespace octopus

#endif
