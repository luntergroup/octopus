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
#include <iterator>
#include <type_traits>

#include "aligned_read.hpp"

namespace Octopus
{
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
    
    template <typename InputIt>
    void transform_reads(InputIt first, InputIt last) const;
    
private:
    std::vector<ReadTransformation> transforms_;
    
    void transform_read(AlignedRead& read) const;
};

// private methods

template <typename InputIt>
void ReadTransform::transform_reads(InputIt first, InputIt last) const
{
    std::for_each(first, last, [this] (auto& read) { transform_read(read); });
}

// non-member methods

namespace detail
{
    template <typename Container>
    void transform_reads(Container& reads, const ReadTransform& transformer, std::true_type)
    {
        transformer.transform_reads(std::begin(reads), std::end(reads));
    }
    
    template <typename ReadMap>
    void transform_reads(ReadMap& reads, const ReadTransform& transformer, std::false_type)
    {
        for (auto& p : reads) {
            transform_reads(p.second, transformer, std::true_type {});
        }
    }
} // namespace detail

template <typename Container>
void transform_reads(Container& reads, const ReadTransform& transformer)
{
    using ValueType = typename std::decay_t<typename Container::value_type>;
    detail::transform_reads(reads, transformer, std::is_same<ValueType, AlignedRead> {});
}

} // namespace Octopus

#endif
