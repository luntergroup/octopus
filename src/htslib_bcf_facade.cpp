//
//  htslib_bcf_facade.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "htslib_bcf_facade.h"

#include <vector>
#include <map>

#include "genomic_region.h"
#include "variant.h"
#include "vcf_record.h"

#include <iostream> // TEST

// public methods

HtslibBcfFacade::HtslibBcfFacade(const fs::path& file_path)
:
the_file_path_ {file_path},
the_file_ {bcf_open(the_file_path_.string().c_str(), "r"), htslib_file_deleter},
the_header_ {(the_file_ != nullptr) ? bcf_hdr_read(the_file_.get()) : nullptr, htslib_bcf_header_deleter}
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
    HtsBcfSrPtr sr {bcf_sr_init(), htslib_bcf_srs_deleter};
    
    bcf_sr_set_regions(sr.get(), to_string(region).c_str(), 0); // must go before bcf_sr_add_reader
    
    if (!bcf_sr_add_reader(sr.get(), the_file_path_.string().c_str())) {
        throw std::runtime_error {"Failed to open file " + the_file_path_.string()};
    }
    
    bcf1_t* record; // points into sr - don't need to delete
    
    int ninfo {};
    int* intinfo {nullptr};
    float* floatinfo {nullptr};
    char** stringinfo {nullptr};
    int* flaginfo {nullptr}; // not actually populated
    
    int nformat {};
    int* intformat {nullptr};
    float* floatformat {nullptr};
    char** stringformat {nullptr};
    
    std::vector<VcfRecord> result {};
    
    while (bcf_sr_next_line(sr.get())) {
        record = bcf_sr_get_line(sr.get(), 0);
        
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
        
        std::map<std::string, std::vector<std::string>> info {};
        
        for (unsigned i {}; i < record->n_info; ++i) {
            const char* key {the_header_->id[BCF_DT_ID][record->d.info[i].key].key};
            auto type = bcf_hdr_id2type(the_header_.get(), BCF_HL_INFO, record->d.info[i].key);
            
            std::vector<std::string> vals {};
            vals.reserve(ninfo);
            
            switch (type) {
                case BCF_HT_INT:
                    if (bcf_get_info_int32(the_header_.get(), record, key, &intinfo, &ninfo) > 0) {
                        std::transform(intinfo, intinfo + ninfo, std::back_inserter(vals), [] (auto v) {
                            return std::to_string(v);
                        });
                    }
                    break;
                case BCF_HT_REAL:
                    if (bcf_get_info_float(the_header_.get(), record, key, &floatinfo, &ninfo) > 0) {
                        std::transform(floatinfo, floatinfo + ninfo, std::back_inserter(vals), [] (auto v) {
                            return std::to_string(v);
                        });
                    }
                    break;
                case BCF_HT_STR:
                    if (bcf_get_info_string(the_header_.get(), record, key, &stringinfo, &ninfo) > 0) {
                        std::for_each(stringinfo, stringinfo + ninfo, [&vals] (const char* str) {
                            vals.emplace_back(str);
                        });
                    }
                    break;
                case BCF_HT_FLAG:
                    vals.emplace_back((bcf_get_info_flag(the_header_.get(), record, key, &flaginfo, &ninfo) == 1) ? "1" : "0");
                    break;
            }
            
            info.emplace(key, std::move(vals));
        }
        
        std::vector<std::string> format {};
        format.reserve(record->n_fmt);
        
        for (unsigned i {}; i < record->n_fmt; ++i) {
            format.emplace_back(the_header_->id[BCF_DT_ID][record->d.fmt[i].id].key);
        }
        
        std::map<std::string, std::vector<std::string>> genotypes {};
        
        format.erase(format.begin() + 1, format.end()); // TEST
        
        for (const auto& key : format) {
            std::cout << key << std::endl;
            
            auto type = bcf_hdr_id2type(the_header_.get(), BCF_HL_FMT, bcf_hdr_id2int(the_header_.get(), BCF_DT_ID, key.c_str()));
            
            std::vector<std::string> vals {};
            vals.reserve(record->n_sample);
            
            switch (type) {
                case BCF_HT_INT:
                    if (bcf_get_format_int32(the_header_.get(), record, key.c_str(), &intformat, &nformat) > 0) {
                        std::transform(intformat, intformat + record->n_sample, std::back_inserter(vals), [] (auto v) {
                            return std::to_string(v);
                        });
                    }
                    break;
                case BCF_HT_REAL:
                    if (bcf_get_format_float(the_header_.get(), record, key.c_str(), &floatformat, &nformat) > 0) {
                        std::transform(floatformat, floatformat + record->n_sample, std::back_inserter(vals), [] (auto v) {
                            return std::to_string(v);
                        });
                    }
                    break;
                case BCF_HT_STR:
                    if (bcf_get_format_string(the_header_.get(), record, key.c_str(), &stringformat, &nformat) > 0) {
                        std::for_each(stringformat, stringformat + record->n_sample, [&vals] (const char* str) {
                            vals.emplace_back(str);
                        });
                    }
                    break;
            }
            
            for (unsigned i {}; i < record->n_sample; ++i) {
                genotypes[the_header_->samples[i]].push_back(std::move(vals[i]));
            }
        }
        
        result.emplace_back(std::move(chrom), pos, std::move(id), std::move(ref), std::move(alt), qual, std::move(filter), std::move(info), std::move(format), std::move(genotypes));
    }
    
    if (intinfo != nullptr)    delete [] intinfo;
    if (floatinfo != nullptr)  delete [] floatinfo;
    if (stringinfo != nullptr) delete [] stringinfo;
    if (flaginfo != nullptr)   delete [] flaginfo;
    
    return result;
}

// private methods
