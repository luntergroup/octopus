// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_writer.hpp"

#include <utility>

#include "basics/aligned_read.hpp"

namespace octopus { namespace io {

ReadWriter::ReadWriter(Path bam_out, Path bam_template)
: path_ {std::move(bam_out)}
, impl_ {std::make_unique<HtslibSamFacade>(path_, std::move(bam_template))}
{}

ReadWriter::ReadWriter(ReadWriter&& other)
{
    std::lock_guard<std::mutex> lock {other.mutex_};
    path_ = std::move(other.path_);
    impl_ = std::move(other.impl_);
}

void swap(ReadWriter& lhs, ReadWriter& rhs) noexcept
{
    if (&lhs == &rhs) return;
    std::lock(lhs.mutex_, rhs.mutex_);
    std::lock_guard<std::mutex> lock_lhs {lhs.mutex_, std::adopt_lock}, lock_rhs {rhs.mutex_, std::adopt_lock};
    using std::swap;
    swap(lhs.path_, rhs.path_);
    swap(lhs.impl_, rhs.impl_);
}

const ReadWriter::Path& ReadWriter::path() const noexcept
{
    return path_;
}

void ReadWriter::write(const AlignedRead& read)
{
    std::lock_guard<std::mutex> lock {mutex_};
    impl_->write(read);
}

ReadWriter& operator<<(ReadWriter& dst, const AlignedRead& read)
{
    dst.write(read);
    return dst;
}

} // namespace io
} // namespace octopus
