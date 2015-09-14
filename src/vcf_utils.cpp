//
//  vcf_utils.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "vcf_utils.h"

unsigned get_field_cardinality(const VcfHeader::KeyType& key, const VcfRecord& record)
{
    return 0;
}

std::vector<VcfType> get_typed_info_values(const VcfHeader& header, const VcfRecord& record,
                              const VcfHeader::KeyType& key)
{
    return get_typed_info_values(header, key, record.get_info_value(key));
}

std::vector<VcfType> get_typed_format_values(const VcfHeader& header, const VcfRecord& record,
                                const VcfRecord::SampleIdType sample, const VcfHeader::KeyType& key)
{
    return get_typed_format_values(header, key, record.get_sample_value(sample, key));
}
