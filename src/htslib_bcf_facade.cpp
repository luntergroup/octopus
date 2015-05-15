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

#include <iostream> // TEST

HtslibBcfFacade::HtslibBcfFacade(const fs::path& the_file_path)
:
the_file_path_ {the_file_path},
the_file_ {bcf_open(the_file_path.string().c_str(), "r"), htslib_file_deleter},
the_header_ {bcf_hdr_read(the_file_.get()), htslib_bcf_header_deleter}
{
    if (the_file_ == nullptr) {
        throw std::runtime_error {"Could not initalise memory for file " + the_file_path_.string()};
    }
    
    if (the_header_ == nullptr) {
        throw std::runtime_error {"Cannot find header for " + the_file_path_.string()};
    }
}

std::vector<Variant> HtslibBcfFacade::fetch_variants(const GenomicRegion& a_region)
{
    HtsBcfSrPtr ptr {bcf_sr_init(), htslib_bcf_srs_deleter};
    
    bcf_sr_set_regions(ptr.get(), to_string(a_region).c_str(), 0);

    if (!bcf_sr_add_reader(ptr.get(), the_file_path_.string().c_str())) {
        throw std::runtime_error {"Failed to open file " + the_file_path_.string()};
    }
    
    std::vector<Variant> result {};
    
    while (bcf_sr_next_line(ptr.get())) {
        bcf1_t* record = bcf_sr_get_line(ptr.get(), 0);
        
        const char* contig = bcf_hdr_id2name(the_header_.get(), record->rid);
        auto begin = static_cast<GenomicRegion::SizeType>(record->pos);
        auto end   = begin + static_cast<GenomicRegion::SizeType>(record->rlen);
        GenomicRegion the_region {contig, begin, end};
        
        Variant v {the_region, record->d.allele[0], record->d.allele[1]};
        
        std::cout << v << std::endl;
    }
    
    return result;
}

void HtslibBcfFacade::set_region(const GenomicRegion& a_region)
{
    
}
