//
//  htslib_bcf_facade.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "htslib_bcf_facade.h"

#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <algorithm> // std::transform, std::for_each
#include <iterator>  // std::cbegin, std::cend, std::next
#include <utility>   // std::move
#include <cstring>   // std::strcpy
#include <cstdint>

#include "genomic_region.h"
#include "vcf_header.h"
#include "vcf_record.h"

#include <iostream> // TEST

char* stringcopy(const std::string& source)
{
    char* result = new char[source.length() + 1];
    std::strcpy(result, source.c_str());
    return result;
}

// public methods

HtslibBcfFacade::HtslibBcfFacade(const fs::path& file_path, const std::string& mode)
:
file_path_ {file_path},
file_ {bcf_open(file_path_.string().c_str(), mode.c_str()), htslib_file_deleter},
header_ {(file_ != nullptr && mode == "r") ? bcf_hdr_read(file_.get()) : bcf_hdr_init(mode.c_str()), htslib_bcf_header_deleter},
samples_ {}
{
    if (mode == "r" && file_ == nullptr) {
        throw std::runtime_error {"Could not initalise memory for file " + file_path_.string()};
    }
    
    if (header_ == nullptr) {
        throw std::runtime_error {"Could not make header for file " + file_path_.string()};
    }
    
    if (mode == "r") {
        for (unsigned sample {}; sample < bcf_hdr_nsamples(header_.get()); ++sample) {
            samples_.emplace_back(header_->samples[sample]);
        }
    }
}

std::unordered_map<std::string, std::string> get_format(bcf_hrec_t* line);

VcfHeader HtslibBcfFacade::fetch_header()
{
    VcfHeader result {};
    
    result.set_file_format(bcf_hdr_get_version(header_.get()));
    result.put_samples(samples_);
    
    for (unsigned i {}; i < header_->nhrec; ++i) {
        const auto record = header_->hrec[i];
        switch (record->type) {
            case BCF_HL_GEN: // key=value
                result.put_field(record->key, record->value);
                break;
            default: // TAG=<A=..,B=..>
                result.put_field(record->key, get_format(record));
                break;
        }
    }
    
    return result;
}

std::vector<VcfRecord> HtslibBcfFacade::fetch_records()
{
    HtsBcfSrPtr sr {bcf_sr_init(), htslib_bcf_srs_deleter};
    
    if (!bcf_sr_add_reader(sr.get(), file_path_.string().c_str())) {
        throw std::runtime_error {"Failed to open file " + file_path_.string()};
    }
    
    return fetch_records(sr);
}

std::vector<VcfRecord> HtslibBcfFacade::fetch_records(const GenomicRegion& region)
{
    HtsBcfSrPtr sr {bcf_sr_init(), htslib_bcf_srs_deleter};
    
    bcf_sr_set_regions(sr.get(), to_string(region).c_str(), 0); // must go before bcf_sr_add_reader
    
    if (!bcf_sr_add_reader(sr.get(), file_path_.string().c_str())) {
        throw std::runtime_error {"Failed to open file " + file_path_.string()};
    }
    
    return fetch_records(sr);
}

auto get_hts_tag_type(const std::string& tag)
{
    const static std::unordered_map<std::string, int> types {
        {"FILTER", BCF_HL_FLT},
        {"INFO", BCF_HL_INFO},
        {"FORMAT", BCF_HL_FMT},
        {"CONTIG", BCF_HL_CTG}
    };
    
    return (types.count(tag) == 1) ? types.at(tag) : BCF_HL_STR;
}

void HtslibBcfFacade::write(const VcfHeader& header)
{
    auto hdr = bcf_hdr_init("w");
    
    bcf_hdr_set_version(hdr, header.get_file_format().c_str());
    
    for (auto key : header.get_fields()) {
        bcf_hrec_t* hrec {new bcf_hrec_t};
        
        hrec->type  = BCF_HL_GEN;
        hrec->key   = stringcopy(key);
        hrec->value = stringcopy(header.get_field(key));
        
        hrec->nkeys = 0;
        hrec->keys  = nullptr;
        hrec->vals  = nullptr;
        
        bcf_hdr_add_hrec(hdr, hrec);
    }
    
    for (auto tag : header.get_tags()) {
        bcf_hrec_t* hrec {new bcf_hrec_t};
        
        hrec->type  = get_hts_tag_type(tag);
        
        hrec->nkeys = 0;
        hrec->keys  = nullptr;
        hrec->vals  = nullptr;
        
        hrec->key   = nullptr;
        hrec->value = nullptr;
        
        bcf_hdr_add_hrec(hdr, hrec);
    }
    
    for (const auto& sample : header.get_samples()) {
        bcf_hdr_add_sample(hdr, sample.c_str());
    }
    
    bcf_hdr_write(file_.get(), hdr);
    
    header_.reset(hdr);
}

void set_chrom(bcf_hdr_t* header, bcf1_t* record, const std::string& chrom);
void set_pos(bcf_hdr_t* header, bcf1_t* record, VcfRecord::SizeType pos);
void set_id(bcf_hdr_t* header, bcf1_t* record, const std::string& id);
void set_alleles(bcf_hdr_t* header, bcf1_t* record, const VcfRecord::SequenceType& ref, const std::vector<VcfRecord::SequenceType>& alts);
void set_qual(bcf_hdr_t* header, bcf1_t* record, VcfRecord::QualityType qual);
void set_filters(bcf_hdr_t* header, bcf1_t* record, const std::vector<std::string>& filters);
void set_info(bcf_hdr_t* header, bcf1_t* record, const VcfRecord& r);

void HtslibBcfFacade::write(const VcfRecord& record)
{
    auto r = bcf_init();
    
    set_chrom(header_.get(), r, record.get_chromosome_name());
    set_pos(header_.get(), r, record.get_position());
    set_id(header_.get(), r, record.get_id());
    set_alleles(header_.get(), r, record.get_ref_allele(), record.get_alt_alleles());
    set_qual(header_.get(), r, record.get_quality());
    set_filters(header_.get(), r, record.get_filters());
    set_info(header_.get(), r, record);
    
    if (record.has_sample_data()) {
        
    }
    
    bcf_write(file_.get(), header_.get(), r);
    
    bcf_destroy(r);
}

// private and non-member methods

std::unordered_map<std::string, std::string> get_format(bcf_hrec_t* line)
{
    std::unordered_map<std::string, std::string> result {};
    result.reserve(line->nkeys);
    
    for (unsigned k {}; k < line->nkeys; ++k) {
        if (std::strcmp(line->keys[k], "IDX") != 0) {
            result.emplace(line->keys[k], line->vals[k]);
        }
    }
    
    return result;
}

auto get_chrom(bcf_hdr_t* header, bcf1_t* record)
{
    return bcf_hdr_id2name(header, record->rid);
}

void set_chrom(bcf_hdr_t* header, bcf1_t* record, const std::string& chrom)
{
    record->rid = bcf_hdr_name2id(header, chrom.c_str());
}

auto get_pos(bcf_hdr_t* header, bcf1_t* record)
{
    return record->pos;
}

void set_pos(bcf_hdr_t* header, bcf1_t* record, VcfRecord::SizeType pos)
{
    record->pos = static_cast<std::uint32_t>(pos);
}

auto get_id(bcf_hdr_t* header, bcf1_t* record)
{
    return record->d.id;
}

void set_id(bcf_hdr_t* header, bcf1_t* record, const std::string& id)
{
    record->d.id = new char[id.length() + 1];
    std::strcpy(record->d.id, id.c_str());
}

auto get_ref(bcf_hdr_t* header, bcf1_t* record)
{
    return record->d.allele[0];
}

void set_alleles(bcf_hdr_t* header, bcf1_t* record, const VcfRecord::SequenceType& ref, const std::vector<VcfRecord::SequenceType>& alts)
{
    char** alleles = new char*[alts.size() + 1];
    
    alleles[0] = new char[ref.length() + 1];
    std::strcpy(alleles[0], ref.c_str());
    
    for (unsigned i {1}; i <= alts.size(); ++i) {
        alleles[i] = new char[alts[i - 1].length() + 1];
        std::strcpy(alleles[i], alts[i - 1].c_str());
    }
    
    bcf_update_alleles(header, record, (const char**) alleles, static_cast<int>(alts.size() + 1));
}

auto get_alt(bcf_hdr_t* header, bcf1_t* record)
{
    auto num_alleles = record->n_allele;
    
    std::vector<VcfRecord::SequenceType> result {};
    result.reserve(num_alleles - 1); // one is the reference
    
    for (unsigned i {1}; i < num_alleles; ++i) {
        result.emplace_back(record->d.allele[i]);
    }
    
    return result;
}

auto get_qual(bcf_hdr_t* header, bcf1_t* record)
{
    return record->qual;
}

void set_qual(bcf_hdr_t* header, bcf1_t* record, VcfRecord::QualityType qual)
{
    record->qual = static_cast<float>(qual);
}

auto get_filter(bcf_hdr_t* header, bcf1_t* record)
{
    std::vector<VcfRecord::KeyType> result {};
    
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

void set_filters(bcf_hdr_t* header, bcf1_t* record, const std::vector<std::string>& filters)
{
    for (const auto& filter : filters) {
        bcf_add_filter(header, record, bcf_hdr_id2int(header, BCF_DT_ID, filter.c_str()));
    }
}

auto get_info(bcf_hdr_t* header, bcf1_t* record)
{
    int ninfo {};
    int* intinfo      {nullptr};
    float* floatinfo  {nullptr};
    char** stringinfo {nullptr};
    int* flaginfo     {nullptr}; // not actually populated
    
    std::unordered_map<VcfRecord::KeyType, std::vector<std::string>> result {};
    result.reserve(record->n_info);
    
    for (unsigned i {}; i < record->n_info; ++i) {
        const char* key {header->id[BCF_DT_ID][record->d.info[i].key].key};
        
        std::vector<std::string> values {};
        values.reserve(ninfo);
        
        switch (bcf_hdr_id2type(header, BCF_HL_INFO, record->d.info[i].key)) {
            case BCF_HT_INT:
                if (bcf_get_info_int32(header, record, key, &intinfo, &ninfo) > 0) {
                    std::transform(intinfo, intinfo + ninfo, std::back_inserter(values), [] (auto v) {
                        return std::to_string(v);
                    });
                }
                break;
            case BCF_HT_REAL:
                if (bcf_get_info_float(header, record, key, &floatinfo, &ninfo) > 0) {
                    std::transform(floatinfo, floatinfo + ninfo, std::back_inserter(values), [] (auto v) {
                        return std::to_string(v);
                    });
                }
                break;
            case BCF_HT_STR:
                if (bcf_get_info_string(header, record, key, &stringinfo, &ninfo) > 0) {
                    std::for_each(stringinfo, stringinfo + ninfo, [&values] (const char* str) {
                        values.emplace_back(str);
                    });
                }
                break;
            case BCF_HT_FLAG:
                values.emplace_back((bcf_get_info_flag(header, record, key, &flaginfo, &ninfo) == 1) ? "1" : "0");
                break;
        }
        
        result.emplace(key, std::move(values));
    }
    
    if (intinfo    != nullptr) delete[] intinfo;
    if (floatinfo  != nullptr) delete[] floatinfo;
    if (stringinfo != nullptr) delete[] stringinfo;
    if (flaginfo   != nullptr) delete[] flaginfo;
    
    return result;
}

void set_info(bcf_hdr_t* header, bcf1_t* record, const VcfRecord& r)
{
    auto keys = r.get_info_keys();
    
    for (const auto& key : keys) {
        auto values     = r.get_info_value(key);
        auto num_values = static_cast<int>(values.size());
        
        switch (bcf_hdr_id2type(header, BCF_HL_INFO, bcf_hdr_id2int(header, BCF_DT_ID, key.c_str()))) {
            case BCF_HT_INT:
            {
                int* vals = new int[num_values];
                std::transform(std::cbegin(values), std::cend(values), vals, [] (const auto& v) {
                    return std::stoi(v);
                });
                bcf_update_info_int32(header, record, key.c_str(), &vals, num_values);
                break;
            }
            case BCF_HT_REAL:
            {
                float* vals = new float[num_values];
                std::transform(std::cbegin(values), std::cend(values), vals, [] (const auto& v) {
                    return std::stof(v);
                });
                bcf_update_info_float(header, record, key.c_str(), &vals, num_values);
                break;
            }
            case BCF_HT_STR:
            {
                char** vals = new char*[num_values];
                std::transform(std::cbegin(values), std::cend(values), vals, [] (const auto& v) {
                    char* s = new char[v.length() + 1];
                    std::strcpy(s, v.c_str());
                    return s;
                });
                bcf_update_info_string(header, record, key.c_str(), (const char**) &vals);
                break;
            }
            case BCF_HT_FLAG:
            {
                bcf_update_info_flag(header, record, key.c_str(), "", values.front() == "1");
                break;
            }
        }
    }
}

bool has_samples(bcf_hdr_t* header)
{
    return bcf_hdr_nsamples(header) > 0;
}

auto get_format(bcf_hdr_t* header, bcf1_t* record)
{
    std::vector<VcfRecord::KeyType> result {};
    result.reserve(record->n_fmt);
    
    for (unsigned i {}; i < record->n_fmt; ++i) {
        result.emplace_back(header->id[BCF_DT_ID][record->d.fmt[i].id].key);
    }
    
    return result;
}

auto get_samples(bcf_hdr_t* header, bcf1_t* record, const std::vector<VcfRecord::KeyType>& format)
{
    int nformat {};
    int* intformat      {nullptr};
    float* floatformat  {nullptr};
    char** stringformat {nullptr};
    
    auto num_samples = record->n_sample;
    
    std::unordered_map<VcfRecord::SampleIdType, std::pair<std::vector<VcfRecord::SequenceType>, bool>> genotypes {};
    genotypes.reserve(num_samples);
    
    if (format.front() == "GT") { // the first key must be GT if present
        int ngt {}, g {};
        int* gt {nullptr};
        
        bcf_get_genotypes(header, record, &gt, &ngt);
        auto ploidy = record->d.fmt->n;
        
        for (unsigned sample {}, i {}; sample < num_samples; ++sample, i += ploidy) {
            std::vector<VcfRecord::SequenceType> alleles {};
            alleles.reserve(ploidy);
            
            for (unsigned p {}; p < ploidy; ++p) {
                g = gt[i + p];
                alleles.emplace_back(bcf_gt_is_missing(g) ? "." : record->d.allele[bcf_gt_allele(g)]);
            }
            
            genotypes.emplace(header->samples[sample], std::make_pair(std::move(alleles), bcf_gt_is_phased(g)));
        }
        
        delete[] gt;
    }
    
    std::unordered_map<VcfRecord::SampleIdType, std::unordered_map<VcfRecord::KeyType, std::vector<std::string>>> other_data {};
    
    auto end = std::cend(format);
    for (auto it = std::next(std::cbegin(format)); it != end; ++it) {
        const auto& key = *it;
        
        std::vector<std::string> values {};
        values.reserve(num_samples);
        
        switch (bcf_hdr_id2type(header, BCF_HL_FMT, bcf_hdr_id2int(header, BCF_DT_ID, key.c_str()))) {
            case BCF_HT_INT:
                if (bcf_get_format_int32(header, record, key.c_str(), &intformat, &nformat) > 0) {
                    std::transform(intformat, intformat + record->n_sample, std::back_inserter(values), [] (auto v) {
                        return std::to_string(v);
                    });
                }
                break;
            case BCF_HT_REAL:
                if (bcf_get_format_float(header, record, key.c_str(), &floatformat, &nformat) > 0) {
                    std::transform(floatformat, floatformat + record->n_sample, std::back_inserter(values), [] (auto v) {
                        return std::to_string(v);
                    });
                }
                break;
            case BCF_HT_STR:
                if (bcf_get_format_string(header, record, key.c_str(), &stringformat, &nformat) > 0) {
                    std::for_each(stringformat, stringformat + record->n_sample, [&values] (const char* str) {
                        values.emplace_back(str);
                    });
                }
                break;
        }
        
        for (unsigned sample {}; sample < num_samples; ++sample) {
            other_data[header->samples[sample]][key].push_back(std::move(values[sample]));
        }
    }
    
    if (intformat    != nullptr) delete[] intformat;
    if (floatformat  != nullptr) delete[] floatformat;
    if (stringformat != nullptr) delete[] stringformat;
    
    return std::make_pair(genotypes, other_data);
}

std::vector<VcfRecord> HtslibBcfFacade::fetch_records(HtsBcfSrPtr& sr)
{
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
