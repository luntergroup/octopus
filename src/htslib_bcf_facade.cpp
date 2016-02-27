//
//  htslib_bcf_facade.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "htslib_bcf_facade.hpp"

#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <utility>
#include <cstring>
#include <cstdint>

#include <boost/filesystem/operations.hpp>

#include "genomic_region.hpp"
#include "vcf_header.hpp"
#include "vcf_record.hpp"

#include <iostream> // TEST

char* convert(const std::string& source)
{
    const auto result = (char*) std::malloc(source.length() + 1);
    source.copy(result, source.length());
    result[source.length()] = '\0';
    return result;
}

std::vector<std::string> get_samples(const bcf_hdr_t* header)
{
    std::vector<std::string> result {};
    
    const auto num_samples = bcf_hdr_nsamples(header);
    
    result.reserve(num_samples);
    
    for (unsigned s {0}; s < num_samples; ++s) {
        result.emplace_back(header->samples[s]);
    }
    
    return result;
}

// public methods

std::string get_hts_mode(const HtslibBcfFacade::Path& file_path, const std::string& mode)
{
    if (!(mode == "r" || mode == "w")) {
        throw std::runtime_error {"invalid mode " + mode + " given to HtslibBcfFacade; must be r or w"};
    }
    
    auto result = "[" + mode + "]";
    
    if (mode == "w" && file_path != "-") { // "-" is for stdout
        auto extension = file_path.extension();
        if (extension == ".bcf") {
            result += "b";
        } else if (extension == ".gz" && file_path.stem().extension() == ".vcf") {
            result += "z";
        }
    }
    
    return result;
}

HtslibBcfFacade::HtslibBcfFacade(Path file_path, const std::string& mode)
:
file_path_ {std::move(file_path)},
file_ {nullptr, HtsFileDeleter {}},
header_ {nullptr, HtsHeaderDeleter {}},
samples_ {}
{
    const auto hts_mode = get_hts_mode(file_path_, mode);
    
    if (mode == "r") {
        if (boost::filesystem::exists(file_path_)) {
            file_.reset(bcf_open(file_path_.c_str(), hts_mode.c_str()));
            
            if (file_ == nullptr) return;
            
            header_.reset(bcf_hdr_read(file_.get()));
            
            if (header_ == nullptr) {
                throw std::runtime_error {"HtslibBcfFacade: could not make header for file " + file_path_.string()};
            }
            
            samples_ = get_samples(header_.get());
        }
    } else {
        file_.reset(bcf_open(file_path_.c_str(), hts_mode.c_str()));
        header_.reset(bcf_hdr_init(hts_mode.c_str()));
    }
}

std::unordered_map<std::string, std::string> get_format(const bcf_hrec_t* line);

bool HtslibBcfFacade::is_header_written() const noexcept
{
    return header_ != nullptr;
}

VcfHeader HtslibBcfFacade::fetch_header() const
{
    VcfHeader::Builder result {};
    
    result.set_file_format(bcf_hdr_get_version(header_.get()));
    result.set_samples(samples_);
    
    std::for_each(header_->hrec, header_->hrec + header_->nhrec,
                  [&result] (const auto record) {
                      switch (record->type) {
                          case BCF_HL_GEN: // key=value
                              result.add_basic_field(record->key, record->value);
                              break;
                          default: // TAG=<A=..,B=..>
                              result.add_structured_field(record->key, get_format(record));
                              break;
                      }
                  });
    
    return result.build_once();
}

size_t HtslibBcfFacade::count_records()
{
    HtsBcfSrPtr sr {bcf_sr_init(), HtsSrsDeleter {}};
    
    if (bcf_sr_add_reader(sr.get(), file_path_.c_str()) != 1) {
        sr.release();
        throw std::runtime_error {"failed to open file " + file_path_.string()};
    }
    
    return count_records(sr);
}

size_t HtslibBcfFacade::count_records(const std::string& contig)
{
    HtsBcfSrPtr sr {bcf_sr_init(), HtsSrsDeleter {}};
    
    if (bcf_sr_set_regions(sr.get(), contig.c_str(), 0) != 0) { // must go before bcf_sr_add_reader
        throw std::runtime_error {"failed to load contig " + contig};
    }
    
    if (bcf_sr_add_reader(sr.get(), file_path_.c_str()) != 1) {
        sr.release();
        throw std::runtime_error {"failed to open file " + file_path_.string()};
    }
    
    return count_records(sr);
}

size_t HtslibBcfFacade::count_records(const GenomicRegion& region)
{
    HtsBcfSrPtr sr {bcf_sr_init(), HtsSrsDeleter {}};
    
    if (bcf_sr_set_regions(sr.get(), to_string(region).c_str(), 0) != 0) { // must go before bcf_sr_add_reader
        throw std::runtime_error {"failed to load region " + to_string(region)};
    }
    
    if (bcf_sr_add_reader(sr.get(), file_path_.c_str()) != 1) {
        sr.release();
        throw std::runtime_error {"failed to open file " + file_path_.string()};
    }
    
    return count_records(sr);
}

std::vector<VcfRecord> HtslibBcfFacade::fetch_records(const Unpack level)
{
    const auto n_records = count_records();
    
    HtsBcfSrPtr sr {bcf_sr_init(), HtsSrsDeleter {}};
    
    if (bcf_sr_add_reader(sr.get(), file_path_.c_str()) != 1) {
        sr.release();
        throw std::runtime_error {"failed to open file " + file_path_.string()};
    }
    
    return fetch_records(sr, level, n_records);
}

std::vector<VcfRecord> HtslibBcfFacade::fetch_records(const std::string& contig, const Unpack level)
{
    const auto n_records = count_records(contig);
    
    HtsBcfSrPtr sr {bcf_sr_init(), HtsSrsDeleter {}};
    
    if (bcf_sr_set_regions(sr.get(), contig.c_str(), 0) != 0) { // must go before bcf_sr_add_reader
        throw std::runtime_error {"failed load contig " + contig};
    }
    
    if (bcf_sr_add_reader(sr.get(), file_path_.c_str()) != 1) {
        sr.release();
        throw std::runtime_error {"failed to open file " + file_path_.string()};
    }
    
    return fetch_records(sr, level, n_records);
}

std::vector<VcfRecord> HtslibBcfFacade::fetch_records(const GenomicRegion& region, const Unpack level)
{
    const auto n_records = count_records(region);
    
    HtsBcfSrPtr sr {bcf_sr_init(), HtsSrsDeleter {}};
    
    const auto region_str = to_string(region);
    
    if (bcf_sr_set_regions(sr.get(), region_str.c_str(), 0) != 0) { // must go before bcf_sr_add_reader
        throw std::runtime_error {"failed load region " + region_str};
    }
    
    if (bcf_sr_add_reader(sr.get(), file_path_.c_str()) != 1) {
        sr.release();
        throw std::runtime_error {"failed to open file " + file_path_.string()};
    }
    
    return fetch_records(sr, level, n_records);
}

auto hts_tag_type(const std::string& tag)
{
    const static std::unordered_map<std::string, int> types {
        {"FILTER", BCF_HL_FLT},
        {"INFO", BCF_HL_INFO},
        {"FORMAT", BCF_HL_FMT},
        {"CONTIG", BCF_HL_CTG}
    };
    
    return (types.count(tag) == 1) ? types.at(tag) : BCF_HL_STR;
}

void HtslibBcfFacade::write_header(const VcfHeader& header)
{
    auto hdr = bcf_hdr_init("w");
    
    bcf_hdr_set_version(hdr, header.get_file_format().c_str());
    
    for (auto& p : header.get_basic_fields()) {
        auto hrec = (bcf_hrec_t*) std::malloc(sizeof(bcf_hrec_t));
        
        hrec->type  = BCF_HL_GEN;
        hrec->key   = convert(p.first);
        hrec->value = convert(p.second);
        
        hrec->nkeys = 0;
        hrec->keys  = nullptr;
        hrec->vals  = nullptr;
        
        bcf_hdr_add_hrec(hdr, hrec);
    }
    
    for (auto& tag : header.get_structured_field_tags()) {
        const auto type = hts_tag_type(tag);
        
        for (auto fields : header.get_structured_fields(tag)) {
            auto hrec = (bcf_hrec_t*) std::malloc(sizeof(bcf_hrec_t));
            
            hrec->type  = type;
            hrec->key   = convert(tag);
            hrec->nkeys = static_cast<int>(fields.size());
            hrec->keys  = (char**) std::malloc(sizeof(char*) * fields.size());
            hrec->vals  = (char**) std::malloc(sizeof(char*) * fields.size());
            
            unsigned i {0};
            
            static const std::vector<std::string> preset_fields {"ID", "Number", "Type", "Source", "Version"};
            
            // explicitly writing these preset fields first to get nice ordering in file
            for (const auto& field : preset_fields) {
                if (fields.count(field) == 1) {
                    hrec->keys[i] = convert(field);
                    hrec->vals[i] = convert(fields[field]);
                    fields.erase(field);
                    ++i;
                }
            }
            
            // the rest of the fields go in whatever order they come in the map
            for (const auto& field : fields) {
                hrec->keys[i] = convert(field.first);
                hrec->vals[i] = convert(field.second);
                ++i;
            }
            
            hrec->value = nullptr;
            
            bcf_hdr_add_hrec(hdr, hrec);
        }
    }
    
    for (const auto& sample : header.get_samples()) {
        bcf_hdr_add_sample(hdr, sample.c_str());
    }
    
    bcf_hdr_write(file_.get(), hdr);
    
    header_.reset(hdr);
    
    samples_ = get_samples(header_.get());
}

void set_chrom(const bcf_hdr_t* header, bcf1_t* record, const std::string& chrom);
void set_pos(bcf1_t* record, VcfRecord::SizeType pos);
void set_id(bcf1_t* record, const std::string& id);
void set_alleles(const bcf_hdr_t* header, bcf1_t* record, const VcfRecord::SequenceType& ref,
                 const std::vector<VcfRecord::SequenceType>& alts);
void set_qual(bcf1_t* record, VcfRecord::QualityType qual);
void set_filters(const bcf_hdr_t* header, bcf1_t* record, const std::vector<std::string>& filters);
void set_info(const bcf_hdr_t* header, bcf1_t* dest, const VcfRecord& source);
void set_samples(const bcf_hdr_t* header, bcf1_t* dest, const VcfRecord& source,
                 const std::vector<std::string>& samples);

void HtslibBcfFacade::write_record(const VcfRecord& record)
{
    const auto& contig = record.get_chromosome_name();
    
    if (bcf_hdr_get_hrec(header_.get(), BCF_HL_CTG, "ID", contig.c_str(), nullptr) == nullptr) {
        throw std::runtime_error {"required contig header line missing for contig " + contig};
    }
    
    auto hts_record = bcf_init();
    
    set_chrom(header_.get(), hts_record, contig);
    set_pos(hts_record, record.get_position());
    set_id(hts_record, record.get_id());
    set_alleles(header_.get(), hts_record, record.get_ref_allele(), record.get_alt_alleles());
    set_qual(hts_record, record.get_quality());
    set_filters(header_.get(), hts_record, record.get_filters());
    set_info(header_.get(), hts_record, record);
    
    if (record.num_samples() > 0) {
        set_samples(header_.get(), hts_record, record, samples_);
    }
    
    bcf_write(file_.get(), header_.get(), hts_record);
    
    bcf_destroy(hts_record);
}

// private and non-member methods

std::unordered_map<std::string, std::string> get_format(const bcf_hrec_t* line)
{
    std::unordered_map<std::string, std::string> result {};
    result.reserve(line->nkeys);
    
    for (unsigned k {0}; k < line->nkeys; ++k) {
        if (std::strcmp(line->keys[k], "IDX") != 0) {
            result.emplace(line->keys[k], line->vals[k]);
        }
    }
    
    return result;
}

auto get_chrom(const bcf_hdr_t* header, const bcf1_t* record)
{
    return bcf_hdr_id2name(header, record->rid);
}

void set_chrom(const bcf_hdr_t* header, bcf1_t* record, const std::string& chrom)
{
    record->rid = bcf_hdr_name2id(header, chrom.c_str());
}

auto get_pos(const bcf1_t* record)
{
    return record->pos;
}

void set_pos(bcf1_t* record, const VcfRecord::SizeType pos)
{
    record->pos = static_cast<std::uint32_t>(pos);
}

auto get_id(const bcf1_t* record)
{
    return record->d.id;
}

void set_id(bcf1_t* record, const std::string& id)
{
    record->d.id = convert(id);
}

auto get_ref(const bcf1_t* record)
{
    return record->d.allele[0];
}

void set_alleles(const bcf_hdr_t* header, bcf1_t* record, const VcfRecord::SequenceType& ref,
                 const std::vector<VcfRecord::SequenceType>& alts)
{
    const auto num_alleles = alts.size() + 1;
    
    std::vector<const char*> alleles(num_alleles);
    
    alleles.front() = ref.c_str();
    
    std::transform(std::begin(alts), std::end(alts), std::next(std::begin(alleles)),
                   [] (const auto& allele) { return allele.c_str(); });
    
    bcf_update_alleles(header, record, alleles.data(), static_cast<int>(num_alleles));
}

auto get_alt(const bcf1_t* record)
{
    const auto num_alleles = record->n_allele;
    
    std::vector<VcfRecord::SequenceType> result {};
    result.reserve(num_alleles - 1); // first is the reference
    
    for (unsigned i {1}; i < num_alleles; ++i) {
        result.emplace_back(record->d.allele[i]);
    }
    
    return result;
}

auto get_qual(const bcf1_t* record)
{
    return record->qual;
}

void set_qual(bcf1_t* record, const VcfRecord::QualityType qual)
{
    record->qual = static_cast<float>(qual);
}

auto get_filter(const bcf_hdr_t* header, const bcf1_t* record)
{
    std::vector<VcfRecord::KeyType> result {};
    result.reserve(record->d.n_flt);
    
    for (unsigned i {0}; i < record->d.n_flt; ++i) {
        result.emplace_back(bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[i]));
    }
    
    return result;
}

void set_filters(const bcf_hdr_t* header, bcf1_t* record, const std::vector<std::string>& filters)
{
    for (const auto& filter : filters) {
        bcf_add_filter(header, record, bcf_hdr_id2int(header, BCF_DT_ID, filter.c_str()));
    }
}

auto get_info(const bcf_hdr_t* header, bcf1_t* record)
{
    int nintinfo {};
    int* intinfo      {nullptr};
    int nfloatinfo {};
    float* floatinfo  {nullptr};
    int nstringinfo {};
    char** stringinfo {nullptr};
    int nflaginfo {};
    int* flaginfo     {nullptr}; // not actually populated
    
    std::unordered_map<VcfRecord::KeyType, std::vector<std::string>> result {};
    result.reserve(record->n_info);
    
    for (unsigned i {0}; i < record->n_info; ++i) {
        const char* key {header->id[BCF_DT_ID][record->d.info[i].key].key};
        
        std::vector<std::string> values {};
        
        switch (bcf_hdr_id2type(header, BCF_HL_INFO, record->d.info[i].key)) {
            case BCF_HT_INT:
                if (bcf_get_info_int32(header, record, key, &intinfo, &nintinfo) > 0) {
                    values.reserve(nintinfo);
                    std::transform(intinfo, intinfo + nintinfo, std::back_inserter(values),
                                   [] (auto v) { return std::to_string(v); });
                }
                break;
            case BCF_HT_REAL:
                if (bcf_get_info_float(header, record, key, &floatinfo, &nfloatinfo) > 0) {
                    values.reserve(nfloatinfo);
                    std::transform(floatinfo, floatinfo + nfloatinfo, std::back_inserter(values),
                                   [] (auto v) { return std::to_string(v); });
                }
                break;
            case BCF_HT_STR:
                if (bcf_get_info_string(header, record, key, &stringinfo, &nstringinfo) > 0) {
                    values.reserve(nstringinfo);
                    std::for_each(stringinfo, stringinfo + nstringinfo,
                                  [&values] (const char* str) { values.emplace_back(str); });
                }
                break;
            case BCF_HT_FLAG:
                values.reserve(1);
                values.emplace_back((bcf_get_info_flag(header, record, key, &flaginfo, &nflaginfo) == 1) ? "1" : "0");
                break;
        }
        
        result.emplace(key, std::move(values));
    }
    
    std::free(intinfo);
    std::free(floatinfo);
    std::free(stringinfo);
    std::free(flaginfo);
    
    return result;
}

void set_info(const bcf_hdr_t* header, bcf1_t* dest, const VcfRecord& source)
{
    for (const auto& key : source.get_info_keys()) {
        const auto& values    = source.get_info_value(key);
        const auto num_values = static_cast<int>(values.size());
        
        switch (bcf_hdr_id2type(header, BCF_HL_INFO, bcf_hdr_id2int(header, BCF_DT_ID, key.c_str()))) {
            case BCF_HT_INT:
            {
                std::vector<int> vals(num_values);
                std::transform(std::cbegin(values), std::cend(values), std::begin(vals),
                               [] (const auto& v) { return std::stoi(v); });
                bcf_update_info_int32(header, dest, key.c_str(), vals.data(), num_values);
                break;
            }
            case BCF_HT_REAL:
            {
                std::vector<float> vals(num_values);
                std::transform(std::cbegin(values), std::cend(values), std::begin(vals),
                               [] (const auto& v) { return std::stof(v); });
                bcf_update_info_float(header, dest, key.c_str(), vals.data(), num_values);
                break;
            }
            case BCF_HT_STR:
            {
                std::vector<const char*> vals(num_values);
                std::transform(std::cbegin(values), std::cend(values), std::begin(vals),
                               [] (const auto& value) { return value.c_str(); });
                bcf_update_info_string(header, dest, key.c_str(), vals.data());
                break;
            }
            case BCF_HT_FLAG:
            {
                bcf_update_info_flag(header, dest, key.c_str(), "", values.empty() || values.front() == "1");
                break;
            }
        }
    }
}

bool has_samples(const bcf_hdr_t* header)
{
    return bcf_hdr_nsamples(header) > 0;
}

auto get_format(const bcf_hdr_t* header, const bcf1_t* record)
{
    std::vector<VcfRecord::KeyType> result {};
    result.reserve(record->n_fmt);
    
    for (unsigned i {0}; i < record->n_fmt; ++i) {
        result.emplace_back(header->id[BCF_DT_ID][record->d.fmt[i].id].key);
    }
    
    return result;
}

auto get_samples(const bcf_hdr_t* header, bcf1_t* record, const std::vector<VcfRecord::KeyType>& format)
{
    const auto num_samples = record->n_sample;
    
    std::unordered_map<VcfRecord::SampleIdType, std::pair<std::vector<VcfRecord::SequenceType>, bool>> genotypes {};
    genotypes.reserve(num_samples);
    
    if (format.front() == "GT") { // the first key must be GT if present
        int ngt {}, g {};
        int* gt {nullptr};
        
        bcf_get_genotypes(header, record, &gt, &ngt);
        const auto ploidy = record->d.fmt->n;
        
        for (unsigned sample {0}, i {}; sample < num_samples; ++sample, i += ploidy) {
            std::vector<VcfRecord::SequenceType> alleles {};
            alleles.reserve(ploidy);
            
            for (unsigned p {0}; p < ploidy; ++p) {
                g = gt[i + p];
                alleles.emplace_back(bcf_gt_is_missing(g) ? "." : record->d.allele[bcf_gt_allele(g)]);
            }
            
            genotypes.emplace(header->samples[sample], std::make_pair(std::move(alleles), bcf_gt_is_phased(g)));
        }
        
        std::free(gt);
    }
    
    std::unordered_map<VcfRecord::SampleIdType, std::unordered_map<VcfRecord::KeyType, std::vector<std::string>>> other_data {};
    
    int nintformat {};
    int* intformat {nullptr};
    int nfloatformat {};
    float* floatformat {nullptr};
    int nstringformat {};
    char** stringformat {nullptr};
    
    for (auto it = std::next(std::cbegin(format)), end = std::cend(format); it != end; ++it) {
        const auto& key = *it;
        
        std::vector<std::string> values {};
        
        switch (bcf_hdr_id2type(header, BCF_HL_FMT, bcf_hdr_id2int(header, BCF_DT_ID, key.c_str()))) {
            case BCF_HT_INT:
                if (bcf_get_format_int32(header, record, key.c_str(), &intformat, &nintformat) > 0) {
                    values.reserve(nintformat);
                    std::transform(intformat, intformat + num_samples, std::back_inserter(values),
                                   [] (auto v) { return std::to_string(v); });
                }
                break;
            case BCF_HT_REAL:
                if (bcf_get_format_float(header, record, key.c_str(), &floatformat, &nfloatformat) > 0) {
                    values.reserve(nfloatformat);
                    std::transform(floatformat, floatformat + num_samples, std::back_inserter(values),
                                   [] (auto v) { return std::to_string(v); });
                }
                break;
            case BCF_HT_STR:
                if (bcf_get_format_string(header, record, key.c_str(), &stringformat, &nstringformat) > 0) {
                    values.reserve(nstringformat);
                    std::for_each(stringformat, stringformat + num_samples,
                                  [&values] (const char* str) { values.emplace_back(str); });
                }
                break;
        }
        
        for (unsigned sample {0}; sample < num_samples; ++sample) {
            other_data[header->samples[sample]][key].push_back(std::move(values[sample]));
        }
    }
    
    std::free(intformat);
    std::free(floatformat);
    std::free(stringformat);
    
    return std::make_pair(std::move(genotypes), std::move(other_data));
}

template <typename T, typename Container>
auto genotype_number(const T& allele, const Container& alleles, const bool is_phased)
{
    if (allele == ".") {
        return (is_phased) ? bcf_gt_missing + 1 : bcf_gt_missing;
    }
    
    const auto it = std::find(std::cbegin(alleles), std::cend(alleles), allele);
    
    const auto allele_num = 2 * static_cast<decltype(bcf_gt_missing)>(std::distance(std::cbegin(alleles), it)) + 2;
    
    return (is_phased) ? allele_num + 1 : allele_num;
}

void set_samples(const bcf_hdr_t* header, bcf1_t* dest, const VcfRecord& source,
                 const std::vector<std::string>& samples)
{
    const auto& alt_alleles = source.get_alt_alleles();
    
    std::vector<VcfRecord::SequenceType> alleles {};
    alleles.reserve(alt_alleles.size() + 1);
    alleles.push_back(source.get_ref_allele());
    alleles.insert(std::end(alleles), std::cbegin(alt_alleles), std::cend(alt_alleles));
    
    const auto num_samples = static_cast<int>(source.num_samples());
    
    if (source.has_genotypes()) {
        const auto ngt = num_samples * static_cast<int>(source.sample_ploidy());
        
        std::vector<int> genotype(ngt);
        
        auto it = std::begin(genotype);
        
        for (const auto& sample : samples) {
            const bool is_phased {source.is_sample_phased(sample)};
            
            const auto& values = source.get_sample_value(sample, "GT");
            
            it = std::transform(std::cbegin(values), std::cend(values), it,
                                [is_phased, &alleles] (const auto& allele) {
                                    return genotype_number(allele, alleles, is_phased);
                                });
        }
        
        bcf_update_genotypes(header, dest, genotype.data(), ngt);
    }
    
    const auto& format = source.get_format();
    
    for (const auto& key : format) {
        if (key == "GT") continue;
        
        const auto num_values = num_samples * static_cast<int>(source.format_cardinality(key));
        
        switch (bcf_hdr_id2type(header, BCF_HL_FMT, bcf_hdr_id2int(header, BCF_DT_ID, key.c_str()))) {
            case BCF_HT_INT:
            {
                std::vector<int> typed_values(num_values);
                auto it = std::begin(typed_values);
                for (const auto& sample : samples) {
                    const auto& values = source.get_sample_value(sample, key);
                    it = std::transform(std::cbegin(values), std::cend(values), it,
                                        [] (const auto& value) { return std::stoi(value); });
                }
                bcf_update_format_int32(header, dest, key.c_str(), typed_values.data(), num_values);
                break;
            }
            case BCF_HT_REAL:
            {
                std::vector<float> typed_values(num_values);
                auto it = std::begin(typed_values);
                for (const auto& sample : samples) {
                    const auto& values = source.get_sample_value(sample, key);
                    it = std::transform(std::cbegin(values), std::cend(values), it,
                                        [] (const auto& value) { return std::stof(value); });
                }
                bcf_update_format_float(header, dest, key.c_str(), typed_values.data(), num_values);
                break;
            }
            case BCF_HT_STR:
            {
                std::vector<const char*> typed_values(num_values);
                auto it = std::begin(typed_values);
                for (const auto& sample : samples) {
                    const auto& values = source.get_sample_value(sample, key);
                    it = std::transform(std::cbegin(values), std::cend(values), it,
                                        [] (const auto& value) { return value.c_str(); });
                }
                bcf_update_format_string(header, dest, key.c_str(), typed_values.data(), num_values);
                break;
            }
        }
    }
}

size_t HtslibBcfFacade::count_records(HtsBcfSrPtr& sr) const
{
    size_t result {0};
    while (bcf_sr_next_line(sr.get())) ++result;
    return result;
}

std::vector<VcfRecord> HtslibBcfFacade::fetch_records(HtsBcfSrPtr& sr, Unpack level,
                                                      const size_t num_records)
{
    std::vector<VcfRecord> result {};
    result.reserve(num_records);
    
    while (bcf_sr_next_line(sr.get())) {
        auto record = bcf_sr_get_line(sr.get(), 0);
        
        bcf_unpack(record, (level == Unpack::All) ? BCF_UN_ALL : BCF_UN_SHR);
        
        auto chrom  = get_chrom(header_.get(), record);
        auto pos    = static_cast<GenomicRegion::SizeType>(get_pos(record));
        auto id     = get_id(record);
        auto ref    = get_ref(record);
        auto alt    = get_alt(record);
        auto qual   = static_cast<VcfRecord::QualityType>(get_qual(record));
        auto filter = get_filter(header_.get(), record);
        auto info   = get_info(header_.get(), record);
        
        if (level == Unpack::All && has_samples(header_.get())) {
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
