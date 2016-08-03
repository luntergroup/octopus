//
//  vcf_header_factory.hpp
//  Octopus
//
//  Created by Daniel Cooke on 05/06/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef vcf_header_factory_hpp
#define vcf_header_factory_hpp

#include <set>
#include <typeindex>
#include <unordered_map>
#include <functional>

#include <io/variant/vcf_header.hpp>

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

#endif /* vcf_header_factory_hpp */
