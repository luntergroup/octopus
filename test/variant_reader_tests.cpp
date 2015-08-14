//
//  variant_reader_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 15/05/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>

#include "test_common.h"
#include "genomic_region.h"
#include "variant.h"
#include "variant_file_reader.h"
#include "variant_file_factory.h"
#include "htslib_bcf_facade.h"
#include "vcf_record.h"

using std::cout;
using std::endl;

BOOST_AUTO_TEST_CASE(can_read_vcf_files)
{
    HtslibBcfFacade vcf_reader {sample_vcf};
    
    GenomicRegion region {"X", 2000000, 2001000};
    
    auto records = vcf_reader.fetch_records(region);
    
    //for (const auto& record : records) cout << record << endl;
}
