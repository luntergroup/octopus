//
//  vcf_utils.h
//  Octopus
//
//  Created by Daniel Cooke on 11/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__vcf_utils__
#define __Octopus__vcf_utils__

#include "vcf_type.h"
#include "vcf_header.h"
#include "vcf_record.h"

unsigned get_field_cardinality(const VcfHeader::KeyType& key, const VcfRecord& record);

std::vector<VcfType> get_typed_info_values(const VcfHeader& header, const VcfRecord& record,
                                           const VcfHeader::KeyType& key);

std::vector<VcfType> get_typed_format_values(const VcfHeader& header, const VcfRecord& record,
                                             const VcfRecord::SampleIdType sample,
                                             const VcfHeader::KeyType& key);

#endif /* defined(__Octopus__vcf_utils__) */
