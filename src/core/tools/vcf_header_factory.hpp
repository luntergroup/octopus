// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef vcf_header_factory_hpp
#define vcf_header_factory_hpp

#include <set>
#include <typeindex>
#include <unordered_map>
#include <functional>

#include "io/variant/vcf_header.hpp"

namespace octopus {

class VcfHeaderFactory
{
public:
    VcfHeaderFactory() = default;
    
    VcfHeaderFactory(const VcfHeaderFactory&)            = default;
    VcfHeaderFactory& operator=(const VcfHeaderFactory&) = default;
    VcfHeaderFactory(VcfHeaderFactory&&)                 = default;
    VcfHeaderFactory& operator=(VcfHeaderFactory&&)      = default;
    
    ~VcfHeaderFactory() = default;
    
    void register_call_type(std::type_index type);
    
    void annotate(VcfHeader::Builder& hb) const;
    
private:
    using Annotator    = std::function<void(VcfHeader::Builder&)>;
    using AnnotatorMap = std::unordered_map<std::type_index, Annotator>;
    
    std::set<std::type_index> call_types_;
    
    static AnnotatorMap annotators_;
};

} // namespace octopus

#endif
