//
//  read_transform.h
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_read_transform_h
#define Octopus_read_transform_h

#include <functional> // std::function
#include <vector>
#include <algorithm>

#include "aligned_read.h"

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
    
    void register_transform(ReadTransformation a_transform);
    template <typename InputIterator>
    void transform_reads(InputIterator first, InputIterator last) const;
    
private:
    std::vector<ReadTransformation> the_transforms_;
    
    void transform_read_(AlignedRead& a_read) const;
};

inline void ReadTransform::register_transform(ReadTransformation a_transform)
{
    the_transforms_.emplace_back(std::move(a_transform));
}

template <typename InputIterator>
void ReadTransform::transform_reads(InputIterator first, InputIterator last) const
{
    std::for_each(first, last, [this] (auto& a_read) { transform_read_(a_read); });
}

inline void ReadTransform::transform_read_(AlignedRead& a_read) const
{
    for (const auto& transform : the_transforms_) {
        transform(a_read);
    }
}

#endif
