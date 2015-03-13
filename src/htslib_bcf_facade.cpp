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

HtslibBcfFacade::HtslibBcfFacade(const fs::path& the_file_path)
:
the_file_path_ {the_file_path},
the_file_ {bcf_sr_init(), htslib_file_deleter},
the_header_ {bcf_hdr_read(the_file_->readers->file), htslib_bcf_header_deleter}
{
    if (the_file_ == nullptr) {
        throw std::runtime_error {"Could not initalise memory for file " + the_file_path_.string()};
    }
    
    if (!bcf_sr_add_reader(the_file_.get(), the_file_path_.string().c_str())) {
        throw std::runtime_error {"Failed to open file " + the_file_path_.string()};
    }
    
    if (the_header_ == nullptr) {
        throw std::runtime_error {"Cannot find header for " + the_file_path_.string()};
    }
}

std::vector<Variant> HtslibBcfFacade::fetch_variants(const GenomicRegion& a_region)
{
    set_region(a_region);
    
    std::vector<Variant> result {};
    
    while (bcf_sr_next_line(the_file_.get())) {
        bcf1_t* record = the_file_->readers[0].buffer[0];
        
        //GenomicRegion the_region {record->rid, }
    }
    
    return result;
}

void HtslibBcfFacade::write_variants(const std::vector<Variant>& some_variants)
{
    // TODO: implement this
}

void HtslibBcfFacade::set_region(const GenomicRegion& a_region)
{
    bcf_sr_set_regions(the_file_.get(), to_string(a_region).c_str(), 0);
}
