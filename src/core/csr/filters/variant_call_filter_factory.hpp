// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef variant_call_filter_factory_hpp
#define variant_call_filter_factory_hpp

#include <memory>

#include "variant_call_filter.hpp"
#include "io/reference/reference_genome.hpp"
#include "readpipe/read_pipe_fwd.hpp"

namespace octopus { namespace csr {

class VariantCallFilterFactory
{
public:
    VariantCallFilterFactory() = default;
    
    VariantCallFilterFactory(const VariantCallFilterFactory&)            = default;
    VariantCallFilterFactory& operator=(const VariantCallFilterFactory&) = default;
    VariantCallFilterFactory(VariantCallFilterFactory&&)                 = default;
    VariantCallFilterFactory& operator=(VariantCallFilterFactory&&)      = default;
    
    ~VariantCallFilterFactory() = default;
    
    std::unique_ptr<VariantCallFilter> make(const ReferenceGenome& reference, const ReadPipe& read_pipe) const;
};

} // namespace csr

using csr::VariantCallFilterFactory;

} // namespace octopus

#endif /* variant_call_filter_factory_hpp */
