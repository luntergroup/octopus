// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "annotated_aligned_read.hpp"

#include <utility>

#include "utils/map_utils.hpp"

namespace octopus {

AnnotatedAlignedRead::AnnotatedAlignedRead(AlignedRead read)
: read_ {std::move(read)}
{}

const GenomicRegion& AnnotatedAlignedRead::mapped_region() const noexcept
{
    return read().mapped_region();
}

const AlignedRead& AnnotatedAlignedRead::read() const noexcept
{
    return read_;
}

AlignedRead& AnnotatedAlignedRead::read() noexcept
{
    return read_;
}

std::vector<AnnotatedAlignedRead::Tag> AnnotatedAlignedRead::tags() const
{
    return extract_keys(annotations_);
}

const AnnotatedAlignedRead::Annotation& AnnotatedAlignedRead::annotation(const Tag& tag) const
{
    return annotations_.at(tag);
}

AnnotatedAlignedRead::Annotation& AnnotatedAlignedRead::annotation(Tag& tag)
{
    return annotations_.at(tag);
}

void AnnotatedAlignedRead::annotate(Tag tag, Annotation annotation)
{
    annotations_.emplace(std::move(tag), std::move(annotation));
}

bool operator==(const AnnotatedAlignedRead& lhs, const AnnotatedAlignedRead& rhs) noexcept
{
    return lhs.read() == rhs.read();
}

bool operator<(const AnnotatedAlignedRead& lhs, const AnnotatedAlignedRead& rhs) noexcept
{
    return lhs.read() < rhs.read();
}

} // namespace octopus
