// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_transformer.hpp"

namespace octopus { namespace readpipe
{

// public methods

void ReadTransformer::register_transform(ReadTransform transform)
{
    transforms_.emplace_back(std::move(transform));
}

unsigned ReadTransformer::num_transforms() const noexcept
{
    return static_cast<unsigned>(transforms_.size());
}

void ReadTransformer::shrink_to_fit() noexcept
{
    transforms_.shrink_to_fit();
}

void ReadTransformer::transform_read(AlignedRead& read) const
{
    for (const auto& transform : transforms_) {
        transform(read);
    }
}

} // namespace readpipe
} // namespace octopus
