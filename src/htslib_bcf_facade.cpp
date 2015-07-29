//
//  htslib_bcf_facade.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "htslib_bcf_facade.h"

#include <vector>

#include "genomic_region.h"
#include "variant.h"
#include "vcf_record.h"

#include <iostream> // TEST

// public methods

HtslibBcfFacade::HtslibBcfFacade(const fs::path& file_path)
:
the_file_path_ {file_path},
the_file_ {bcf_open(the_file_path_.string().c_str(), "r"), htslib_file_deleter},
the_header_ {bcf_hdr_read(the_file_.get()), htslib_bcf_header_deleter}
{
    if (the_file_ == nullptr) {
        throw std::runtime_error {"Could not initalise memory for file " + the_file_path_.string()};
    }
    
    if (the_header_ == nullptr) {
        throw std::runtime_error {"Cannot find header for " + the_file_path_.string()};
    }
}

std::vector<Variant> HtslibBcfFacade::fetch_variants(const GenomicRegion& region)
{
    HtsBcfSrPtr ptr {bcf_sr_init(), htslib_bcf_srs_deleter};
    
    bcf_sr_set_regions(ptr.get(), to_string(region).c_str(), 0);

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

std::vector<VcfRecord> HtslibBcfFacade::fetch_records(const GenomicRegion& region)
{
    HtsBcfSrPtr ptr {bcf_sr_init(), htslib_bcf_srs_deleter};
    
    bcf_sr_set_regions(ptr.get(), to_string(region).c_str(), 0);
    
    if (!bcf_sr_add_reader(ptr.get(), the_file_path_.string().c_str())) {
        throw std::runtime_error {"Failed to open file " + the_file_path_.string()};
    }
    
    std::vector<VcfRecord> result {};
    
    while (bcf_sr_next_line(ptr.get())) {
        bcf1_t* record = bcf_sr_get_line(ptr.get(), 0);
        
        bcf_unpack(record, BCF_UN_ALL);
        
        const char* chrom = bcf_hdr_id2name(the_header_.get(), record->rid);
        auto pos = static_cast<GenomicRegion::SizeType>(record->pos);
        std::string id {record->d.id};
        std::string ref {record->d.allele[0]};
        
        std::vector<std::string> alt {};
        alt.reserve(record->n_allele - 1);
        for (unsigned i {1}; i < record->n_allele; ++i) {
            alt.emplace_back(record->d.allele[i]);
        }
        
        auto qual = static_cast<VcfRecord::QualityType>(record->qual);
        
        std::vector<std::string> filter {};
        if (record->d.n_flt == 0) {
            filter.reserve(1);
            filter.emplace_back("PASS");
        } else {
            filter.reserve(record->d.n_flt);
            
            for (unsigned i {}; i < record->d.n_flt; ++i) {
                filter.emplace_back(bcf_hdr_int2id(the_header_.get(), BCF_DT_ID, record->d.flt[i]));
            }
        }
        
        std::vector<std::string> info {};
        info.reserve(record->n_info);
        
        for (unsigned i {}; i < record->n_info; ++i) {
            std::string key {the_header_->id[BCF_DT_ID][record->d.info[i].key].key};
            char** p; int* in;
            int r = bcf_get_info_string(the_header_.get(), record, key.c_str(), p, in);
            std::cout << r << " " << (p == nullptr) << " " << (in == nullptr) << std::endl;
        }
        
        result.emplace_back(std::move(chrom), pos, std::move(id), std::move(ref), std::move(alt), qual);
    }
    
    return result;
}

// private methods

void HtslibBcfFacade::set_region(const GenomicRegion& a_region)
{
    
}
