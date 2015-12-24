//
//  read_transform.hpp
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_read_transform_hpp
#define Octopus_read_transform_hpp

#include <vector>
#include <functional>
#include <algorithm>

#include "aligned_read.hpp"

namespace Octopus {

class ReadTransform
{
public:
    using ReadTransformation = std::function<void(AlignedRead&)>;
    
    ReadTransform()  = default;
    ~ReadTransform() = default;
    
    ReadTransform(const ReadTransform&)            = default;
    ReadTransform& operator=(const ReadTransform&) = default;
    ReadTransform(ReadTransform&&)                 = default;
    ReadTransform& operator=(ReadTransform&&)      = default;
    
    void register_transform(ReadTransformation transform);
    
    unsigned num_transforms() const noexcept;
    
    template <typename InputIterator>
    void transform_reads(InputIterator first, InputIterator last) const;
    
private:
    std::vector<ReadTransformation> transforms_;
    
    void transform_read_(AlignedRead& read) const;
};

// public methods

inline void ReadTransform::register_transform(ReadTransformation transform)
{
    transforms_.emplace_back(std::move(transform));
}

inline unsigned ReadTransform::num_transforms() const noexcept
{
    return static_cast<unsigned>(transforms_.size());
}

inline void ReadTransform::transform_read_(AlignedRead& read) const
{
    for (const auto& transform : transforms_) {
        transform(read);
    }
}

// private methods

template <typename InputIterator>
void ReadTransform::transform_reads(InputIterator first, InputIterator last) const
{
    std::for_each(first, last, [this] (auto& read) { transform_read_(read); });
}

} // namespace Octopus

#endif
