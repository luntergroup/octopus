//
//  read_transformer.cpp
//  octopus
//
//  Created by Daniel Cooke on 10/01/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

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
