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
#include <stdexcept>
#include <algorithm> // std::transform, std::for_each
#include <iterator>  // std::cbegin, std::cend, std::next
#include <utility>   // std::move

#include "genomic_region.h"
#include "variant.h"
#include "vcf_record.h"

#include <iostream> // TEST

// public methods

HtslibBcfFacade::HtslibBcfFacade(const fs::path& file_path)
:
file_path_ {file_path},
file_ {bcf_open(file_path_.string().c_str(), "r"), htslib_file_deleter},
header_ {(file_ != nullptr) ? bcf_hdr_read(file_.get()) : nullptr, htslib_bcf_header_deleter}
{
    if (file_ == nullptr) {
        throw std::runtime_error {"Could not initalise memory for file " + file_path_.string()};
    }
    
    if (header_ == nullptr) {
        throw std::runtime_error {"Cannot find header for " + file_path_.string()};
    }
}

auto get_chrom(bcf_hdr_t* header, bcf1_t* record)
{
    return bcf_hdr_id2name(header, record->rid);
}

auto get_pos(bcf_hdr_t* header, bcf1_t* record)
{
    return record->pos;
}

auto get_id(bcf_hdr_t* header, bcf1_t* record)
{
    return record->d.id;
}

auto get_ref(bcf_hdr_t* header, bcf1_t* record)
{
    return record->d.allele[0];
}

auto get_alt(bcf_hdr_t* header, bcf1_t* record)
{
    std::vector<std::string> result {};
    result.reserve(record->n_allele - 1);
    
    for (unsigned i {1}; i < record->n_allele; ++i) {
        result.emplace_back(record->d.allele[i]);
    }
    
    return result;
}

auto get_qual(bcf_hdr_t* header, bcf1_t* record)
{
    return record->qual;
}

auto get_filter(bcf_hdr_t* header, bcf1_t* record)
{
    std::vector<std::string> result {};
    
    if (record->d.n_flt == 0) {
        result.reserve(1);
        result.emplace_back("PASS");
    } else {
        result.reserve(record->d.n_flt);
        for (unsigned i {}; i < record->d.n_flt; ++i) {
            result.emplace_back(bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[i]));
        }
    }
    
    return result;
}

auto get_info(bcf_hdr_t* header, bcf1_t* record)
{
    int ninfo {};
    int* intinfo {nullptr};
    float* floatinfo {nullptr};
    char** stringinfo {nullptr};
    int* flaginfo {nullptr}; // not actually populated
    
    std::map<std::string, std::vector<std::string>> result {};
    
    for (unsigned i {}; i < record->n_info; ++i) {
        const char* key {header->id[BCF_DT_ID][record->d.info[i].key].key};
        auto type = bcf_hdr_id2type(header, BCF_HL_INFO, record->d.info[i].key);
        
        std::vector<std::string> vals {};
        vals.reserve(ninfo);
        
        switch (type) {
            case BCF_HT_INT:
                if (bcf_get_info_int32(header, record, key, &intinfo, &ninfo) > 0) {
                    std::transform(intinfo, intinfo + ninfo, std::back_inserter(vals), [] (auto v) {
                        return std::to_string(v);
                    });
                }
                break;
            case BCF_HT_REAL:
                if (bcf_get_info_float(header, record, key, &floatinfo, &ninfo) > 0) {
                    std::transform(floatinfo, floatinfo + ninfo, std::back_inserter(vals), [] (auto v) {
                        return std::to_string(v);
                    });
                }
                break;
            case BCF_HT_STR:
                if (bcf_get_info_string(header, record, key, &stringinfo, &ninfo) > 0) {
                    std::for_each(stringinfo, stringinfo + ninfo, [&vals] (const char* str) {
                        vals.emplace_back(str);
                    });
                }
                break;
            case BCF_HT_FLAG:
                vals.emplace_back((bcf_get_info_flag(header, record, key, &flaginfo, &ninfo) == 1) ? "1" : "0");
                break;
        }
        
        result.emplace(key, std::move(vals));
    }
    
    if (intinfo    != nullptr) delete [] intinfo;
    if (floatinfo  != nullptr) delete [] floatinfo;
    if (stringinfo != nullptr) delete [] stringinfo;
    if (flaginfo   != nullptr) delete [] flaginfo;
    
    return result;
}

bool has_samples(bcf_hdr_t* header)
{
    return bcf_hdr_nsamples(header) > 0;
}

auto get_format(bcf_hdr_t* header, bcf1_t* record)
{
    std::vector<std::string> result {};
    result.reserve(record->n_fmt);
    
    for (unsigned i {}; i < record->n_fmt; ++i) {
        result.emplace_back(header->id[BCF_DT_ID][record->d.fmt[i].id].key);
    }
    
    return result;
}

auto get_samples(bcf_hdr_t* header, bcf1_t* record, const std::vector<std::string>& format)
{
    int nformat {};
    int* intformat {nullptr};
    float* floatformat {nullptr};
    char** stringformat {nullptr};
    
    std::map<std::string, std::pair<std::vector<std::string>, bool>> genotypes {};
    
    if (format.front() == "GT") { // the first key must be GT if present
        int ngt {};
        int *gt {nullptr};
        
        bcf_get_genotypes(header, record, &gt, &ngt);
        auto ploidy = record->d.fmt->n;
        
        for (int i {}, s {}; i < ngt; i += ploidy, ++s) {
            std::vector<std::string> alleles {};
            alleles.reserve(ploidy);
            
            for (int p = 0; p < ploidy; ++p) {
                if (bcf_gt_is_missing(gt[i + p])) {
                    alleles.emplace_back(".");
                } else {
                    alleles.emplace_back(record->d.allele[bcf_gt_allele(gt[i + p])]);
                }
            }
            
            auto is_phased = static_cast<bool>(bcf_gt_is_phased(gt[i + ploidy - 1]));
            
            genotypes[header->samples[s]] = std::make_pair(std::move(alleles), is_phased);
        }
        
        delete [] gt;
    }
    
    std::map<std::string, std::map<std::string, std::vector<std::string>>> other_data {};
    
    auto end = std::cend(format);
    for (auto it = std::next(std::cbegin(format)); it != end; ++it) {
        const auto& key = *it;
        
        std::vector<std::string> vals {};
        vals.reserve(record->n_sample);
        
        auto type = bcf_hdr_id2type(header, BCF_HL_FMT, bcf_hdr_id2int(header, BCF_DT_ID, key.c_str()));
        
        switch (type) {
            case BCF_HT_INT:
                if (bcf_get_format_int32(header, record, key.c_str(), &intformat, &nformat) > 0) {
                    std::transform(intformat, intformat + record->n_sample, std::back_inserter(vals), [] (auto v) {
                        return std::to_string(v);
                    });
                }
                break;
            case BCF_HT_REAL:
                if (bcf_get_format_float(header, record, key.c_str(), &floatformat, &nformat) > 0) {
                    std::transform(floatformat, floatformat + record->n_sample, std::back_inserter(vals), [] (auto v) {
                        return std::to_string(v);
                    });
                }
                break;
            case BCF_HT_STR:
                if (bcf_get_format_string(header, record, key.c_str(), &stringformat, &nformat) > 0) {
                    std::for_each(stringformat, stringformat + record->n_sample, [&vals] (const char* str) {
                        vals.emplace_back(str);
                    });
                }
                break;
        }
        
        for (unsigned s {}; s < record->n_sample; ++s) {
            other_data[header->samples[s]][key].push_back(std::move(vals[s]));
        }
    }
    
    if (intformat    != nullptr) delete [] intformat;
    if (floatformat  != nullptr) delete [] floatformat;
    if (stringformat != nullptr) delete [] stringformat;
    
    return std::make_pair(genotypes, other_data);
}

std::vector<VcfRecord> HtslibBcfFacade::fetch_records(const GenomicRegion& region)
{
    HtsBcfSrPtr sr {bcf_sr_init(), htslib_bcf_srs_deleter};
    
    bcf_sr_set_regions(sr.get(), to_string(region).c_str(), 0); // must go before bcf_sr_add_reader
    
    if (!bcf_sr_add_reader(sr.get(), file_path_.string().c_str())) {
        throw std::runtime_error {"Failed to open file " + file_path_.string()};
    }
    
    bcf1_t* record; // points into sr - don't need to delete
    
    std::vector<VcfRecord> result {};
    
    while (bcf_sr_next_line(sr.get())) {
        record = bcf_sr_get_line(sr.get(), 0);
        
        bcf_unpack(record, BCF_UN_ALL);
        
        auto chrom  = get_chrom(header_.get(), record);
        auto pos    = static_cast<GenomicRegion::SizeType>(get_pos(header_.get(), record));
        auto id     = get_id(header_.get(), record);
        auto ref    = get_ref(header_.get(), record);
        auto alt    = get_alt(header_.get(), record);
        auto qual   = static_cast<VcfRecord::QualityType>(get_qual(header_.get(), record));
        auto filter = get_filter(header_.get(), record);
        auto info   = get_info(header_.get(), record);
        
        if (has_samples(header_.get())) {
            auto format  = get_format(header_.get(), record);
            auto samples = get_samples(header_.get(), record, format);
            result.emplace_back(chrom, pos, std::move(id), std::move(ref), std::move(alt), qual,
                                std::move(filter), std::move(info), std::move(format),
                                std::move(samples.first), std::move(samples.second));
        } else {
            result.emplace_back(chrom, pos, std::move(id), std::move(ref), std::move(alt), qual,
                                std::move(filter), std::move(info));
        }
    }
    
    return result;
}

// private methods
