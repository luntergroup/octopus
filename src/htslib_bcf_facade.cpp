//
//  htslib_bcf_facade.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "htslib_bcf_facade.h"

#include "genomic_region.h"
#include "variant.h"

HtslibBcfFacade::HtslibBcfFacade(const fs::path& variant_file_path)
: the_bcf_file_path_ {variant_file_path},
the_file_ {hts_open(variant_file_path.string().c_str(), "r"), htslib_file_deleter},
the_header_ {bcf_hdr_read(the_file_.get()), htslib_bcf_header_deleter}
{}

std::vector<Variant> HtslibBcfFacade::fetch_variants(const GenomicRegion& a_region)
{
    std::vector<Variant> result {};
    return result;
}

void HtslibBcfFacade::write_variants(const std::vector<Variant>& some_variants)
{
    // TODO: implement this
}
