// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef annotated_aligned_read_hpp
#define annotated_aligned_read_hpp

#include <string>
#include <vector>
#include <unordered_map>

#include "concepts/mappable.hpp"
#include "concepts/comparable.hpp"
#include "basics/genomic_region.hpp"
#include "basics/aligned_read.hpp"
#include "utils/memory_footprint.hpp"

namespace octopus {

class AnnotatedAlignedRead : public Mappable<AnnotatedAlignedRead>
{
public:
    using Tag = std::string;
    using Annotation = std::string;
    
    AnnotatedAlignedRead() = delete;
    
    AnnotatedAlignedRead(AlignedRead read);
    
    AnnotatedAlignedRead(const AnnotatedAlignedRead&)            = default;
    AnnotatedAlignedRead& operator=(const AnnotatedAlignedRead&) = default;
    AnnotatedAlignedRead(AnnotatedAlignedRead&&)                 = default;
    AnnotatedAlignedRead& operator=(AnnotatedAlignedRead&&)      = default;
    
    ~AnnotatedAlignedRead() = default;
    
    const GenomicRegion& mapped_region() const noexcept;
    
    const AlignedRead& read() const noexcept;
    AlignedRead& read() noexcept;
    
    std::vector<Tag> tags() const;
    
    const Annotation& annotation(const Tag& tag) const;
    Annotation& annotation(Tag& tag);
    
    void annotate(Tag tag, Annotation annotation);

private:
    using AnnotationMap = std::unordered_map<Tag, Annotation>;
    
    AlignedRead read_;
    AnnotationMap annotations_;
};

MemoryFootprint footprint(const AnnotatedAlignedRead& read) noexcept;

bool operator==(const AnnotatedAlignedRead& lhs, const AnnotatedAlignedRead& rhs) noexcept;
bool operator<(const AnnotatedAlignedRead& lhs, const AnnotatedAlignedRead& rhs) noexcept;

} // namespace octopus

#endif
