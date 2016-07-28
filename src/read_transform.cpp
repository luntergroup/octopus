//
//  read_transform.cpp
//  Octopus
//
//  Created by Daniel Cooke on 10/01/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "read_transform.hpp"

namespace octopus {
    // public methods
    
    void ReadTransform::register_transform(ReadTransformation transform)
    {
        transforms_.emplace_back(std::move(transform));
    }
    
    unsigned ReadTransform::num_transforms() const noexcept
    {
        return static_cast<unsigned>(transforms_.size());
    }
    
    void ReadTransform::shrink_to_fit() noexcept
    {
        transforms_.shrink_to_fit();
    }
    
    void ReadTransform::transform_read(AlignedRead& read) const
    {
        for (const auto& transform : transforms_) {
            transform(read);
        }
    }
} // namespace octopus
