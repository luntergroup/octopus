//
//  vcf_header_factory.cpp
//  Octopus
//
//  Created by Daniel Cooke on 05/06/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "vcf_header_factory.hpp"

#include <typeinfo>

#include "call_types.hpp"

namespace octopus
{
    VcfHeaderFactory::AnnotatorMap VcfHeaderFactory::annotators_ =
    {
        {std::type_index(typeid(GermlineVariantCall)), [] (auto& hb) {
            hb.add_info("DMP", "1", "Float", "Dummy model posterior");
        }},
        {std::type_index(typeid(ReferenceCall)), [] (auto& hb) {
            hb.add_info("DMP", "1", "Float", "Dummy model posterior");
        }},
        {std::type_index(typeid(SomaticCall)), [] (auto& hb) {
            hb.add_format("SCR", "2", "Float", "99% credible region of the somatic allele frequency");
            hb.add_info("DMP", "1", "Float", "Dummy model posterior");
        }}
    };
    
    void VcfHeaderFactory::register_call_type(std::type_index type)
    {
        call_types_.insert(type);
    }
    
    void VcfHeaderFactory::annotate(VcfHeader::Builder &hb) const
    {
        for (const auto& type : call_types_) {
            annotators_.at(type)(hb);
        }
    }
} // namespace octopus
