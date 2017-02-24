// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_transformer_hpp
#define read_transformer_hpp

#include <vector>
#include <functional>
#include <algorithm>
#include <iterator>
#include <type_traits>

#include "basics/aligned_read.hpp"
#include "containers/mappable_flat_multi_set.hpp"

namespace octopus { namespace readpipe {

class ReadTransformer
{
    using ReadReferenceVector = std::vector<std::reference_wrapper<AlignedRead>>;
    
public:
    using ReadTransform     = std::function<void(AlignedRead&)>;
    using ReadTemplate      = ReadReferenceVector;
    using TemplateTransform = std::function<void(ReadTemplate&)>;
    
    ReadTransformer() = default;
    
    ReadTransformer(const ReadTransformer&)            = default;
    ReadTransformer& operator=(const ReadTransformer&) = default;
    ReadTransformer(ReadTransformer&&)                 = default;
    ReadTransformer& operator=(ReadTransformer&&)      = default;
    
    ~ReadTransformer() = default;
    
    void add(ReadTransform transform);
    void add(TemplateTransform transform);
    
    unsigned num_transforms() const noexcept;
    
    void shrink_to_fit() noexcept;
    
    template <typename ForwardIt>
    void transform_reads(ForwardIt first, ForwardIt last) const;
        
private:
    std::vector<ReadTransform> read_transforms_;
    std::vector<TemplateTransform> template_transforms_;
    
    void transform_read(AlignedRead& read) const;
    template <typename ForwardIt>
    auto make_references(ForwardIt first, ForwardIt last) const;
    void transform(ReadReferenceVector& reads) const;
    void transform_templates(ReadReferenceVector& reads) const;
    void transform_templates(ReadReferenceVector::iterator first, ReadReferenceVector::iterator last) const;
    void transform_template(ReadTemplate& read_template) const;
};

// private methods

template <typename ForwardIt>
auto ReadTransformer::make_references(ForwardIt first, ForwardIt last) const
{
    ReadReferenceVector result {};
    result.reserve(std::distance(first, last));
    std::copy(first, last, std::back_inserter(result));
    return result;
}

template <typename ForwardIt>
void ReadTransformer::transform_reads(ForwardIt first, ForwardIt last) const
{
    if (!read_transforms_.empty()) {
        std::for_each(first, last, [this] (AlignedRead& read) { transform_read(read); });
    }
    if (!template_transforms_.empty()) {
        auto read_references = make_references(first, last);
        transform(read_references);
    }
}

// non-member methods

namespace detail {

template <typename Container>
void transform_reads(Container& reads, const ReadTransformer& transformer, std::true_type)
{
    transformer.transform_reads(std::begin(reads), std::end(reads));
}

template <typename ReadMap>
void transform_reads(ReadMap& reads, const ReadTransformer& transformer, std::false_type)
{
    for (auto& p : reads) {
        transform_reads(p.second, transformer, std::true_type {});
    }
}

} // namespace detail

template <typename Container>
void transform_reads(Container& reads, const ReadTransformer& transformer)
{
    detail::transform_reads(reads, transformer, std::is_same<typename Container::value_type, AlignedRead> {});
}

} // namespace readpipe
} // namespace octopus

#endif
