// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_transformer.hpp"

#include <iostream>

namespace octopus { namespace readpipe {

// public methods

void ReadTransformer::add(ReadTransform transform)
{
    read_transforms_.push_back(std::move(transform));
}

void ReadTransformer::add(TemplateTransform transform)
{
    template_transforms_.push_back(std::move(transform));
}

unsigned ReadTransformer::num_transforms() const noexcept
{
    return static_cast<unsigned>(read_transforms_.size() + template_transforms_.size());
}

void ReadTransformer::shrink_to_fit() noexcept
{
    read_transforms_.shrink_to_fit();
    template_transforms_.shrink_to_fit();
}

void ReadTransformer::transform_read(AlignedRead& read) const
{
    for (const auto& transform : read_transforms_) {
        transform(read);
    }
}

struct ReadTemplateSorter
{
    bool operator()(const AlignedRead& lhs, const AlignedRead& rhs) const noexcept
    {
        return lhs.name() < rhs.name();
    }
};

template <typename Container>
void group_templates(Container& reads)
{
    std::sort(std::begin(reads), std::end(reads), ReadTemplateSorter {});
}

void ReadTransformer::transform_template(ReadTemplate& read_template) const
{
    for (const auto& transform : template_transforms_) {
        transform(read_template);
    }
}

struct ReadTemplateEqual
{
    ReadTemplateEqual(const AlignedRead& read) : read_ {read} {}
    
    bool operator()(const AlignedRead& other) const noexcept
    {
        return other.name() == read_.name();
    }

private:
    const AlignedRead& read_;
};

void ReadTransformer::transform_templates(ReadReferenceVector::iterator first, ReadReferenceVector::iterator last) const
{
    ReadReferenceVector read_template {};
    read_template.reserve(3);
    for (; first != last;) {
        const auto next = std::find_if_not(std::next(first), last, ReadTemplateEqual {*first});
        read_template.assign(first, next);
        transform_template(read_template);
        first = next;
    }
}

void ReadTransformer::transform_templates(ReadReferenceVector& reads) const
{
    transform_templates(std::begin(reads), std::end(reads));
}

void ReadTransformer::transform(ReadReferenceVector& reads) const
{
    group_templates(reads);
    transform_templates(reads);
}

} // namespace readpipe
} // namespace octopus
