// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "aligned_template.hpp"

namespace octopus {

AlignedTemplate::AlignedTemplate(const AlignedRead& read)
: template_ {read}
, region_ {read.mapped_region()}
{}

AlignedTemplate::AlignedTemplate(std::initializer_list<MappableReferenceWrapper<const AlignedRead>> reads)
: template_ {reads}
{
    if (template_.empty()) {
        throw;
    } else {
        region_ = encompassing_region(template_);
    }
}

const GenomicRegion& AlignedTemplate::mapped_region() const noexcept
{
    return region_;
}

std::size_t AlignedTemplate::size() const noexcept
{
    return template_.size();
}

const AlignedRead& AlignedTemplate::operator[](std::size_t idx) const noexcept
{
    return template_[idx].get();
}

bool operator==(const AlignedTemplate& lhs, const AlignedTemplate& rhs) noexcept
{
    return lhs.template_ == rhs.template_;
}

bool operator<(const AlignedTemplate& lhs, const AlignedTemplate& rhs) noexcept
{
    return lhs.template_ < rhs.template_;
}

bool is_rightmost_segment(const AlignedRead& read) noexcept
{
    return !read.is_marked_multiple_segment_template() || (is_reverse_strand(read) && read.is_marked_last_template_segment());
}

} // namespace octopus
