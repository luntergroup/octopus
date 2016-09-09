// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "htslib_bcf_facade.hpp"

#include <vector>
#include <array>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>
#include <utility>
#include <cstdlib>
#include <cstring>
#include <cstdint>

#include <boost/filesystem/operations.hpp>
#include <boost/optional.hpp>
#include <boost/container/small_vector.hpp>

#include "basics/genomic_region.hpp"
#include "utils/string_utils.hpp"
#include "vcf_spec.hpp"
#include "vcf_header.hpp"
#include "vcf_record.hpp"

#include <iostream> // TEST
#include <vcf.h>

namespace octopus {

namespace {
    static const std::string vcfMissingValue {vcfspec::missingValue};
    
    namespace bc = boost::container;
}

char* convert(const std::string& source)
{
    const auto result = (char*) std::malloc(source.length() + 1);
    source.copy(result, source.length());
    result[source.length()] = '\0';
    return result;
}

std::vector<std::string> extract_samples(const bcf_hdr_t* header)
{
    std::vector<std::string> result {};
    const auto num_samples = static_cast<unsigned>(bcf_hdr_nsamples(header));
    result.reserve(num_samples);
    for (unsigned s {0}; s < num_samples; ++s) {
        result.emplace_back(header->samples[s]);
    }
    return result;
}

// public methods

std::string get_hts_mode(const HtslibBcfFacade::Path& file_path, HtslibBcfFacade::Mode mode)
{
    std::string result {"["};
    using Mode = HtslibBcfFacade::Mode;
    if (mode == Mode::read) {
        result += "r]";
    } else {
        result += 'w';
        auto extension = file_path.extension();
        if (extension == ".bcf") {
            result += "b";
        } else if (extension == ".gz" && file_path.stem().extension() == ".vcf") {
            result += "z";
        }
    }
    return result;
}

HtslibBcfFacade::HtslibBcfFacade()
: file_path_ {}
, file_ {bcf_open("-", "[w]"), HtsFileDeleter {}}
, header_ {bcf_hdr_init("[w]"), HtsHeaderDeleter {}}
, samples_ {}
{}

HtslibBcfFacade::HtslibBcfFacade(Path file_path, Mode mode)
: file_path_ {std::move(file_path)}
, file_ {nullptr, HtsFileDeleter {}}
, header_ {nullptr, HtsHeaderDeleter {}}
, samples_ {}
{
    const auto hts_mode = get_hts_mode(file_path_, mode);
    if (mode == Mode::read) {
        if (boost::filesystem::exists(file_path_)) {
            file_.reset(bcf_open(file_path_.c_str(), hts_mode.c_str()));
            if (file_ == nullptr) return;
            header_.reset(bcf_hdr_read(file_.get()));
            if (header_ == nullptr) {
                throw std::runtime_error {"HtslibBcfFacade: could not make header for file " + file_path_.string()};
            }
            samples_ = extract_samples(header_.get());
        } else {
            throw std::runtime_error {"HtslibBcfFacade: " + file_path_.string() + " does not exist"};
        }
    } else {
        file_.reset(bcf_open(file_path_.c_str(), hts_mode.c_str()));
        header_.reset(bcf_hdr_init(hts_mode.c_str()));
    }
}

std::unordered_map<std::string, std::string> extract_format(const bcf_hrec_t* line);

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
                              result.add_structured_field(record->key, extract_format(record));
                              break;
                      }
                  });
    
    return result.build_once();
}

std::size_t HtslibBcfFacade::count_records() const
{
    HtsBcfSrPtr sr {bcf_sr_init(), HtsSrsDeleter {}};
    
    if (bcf_sr_add_reader(sr.get(), file_path_.c_str()) != 1) {
        sr.release();
        throw std::runtime_error {"failed to open file " + file_path_.string()};
    }
    
    return count_records(sr);
}

std::size_t HtslibBcfFacade::count_records(const std::string& contig) const
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

std::size_t HtslibBcfFacade::count_records(const GenomicRegion& region) const
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

HtslibBcfFacade::RecordIteratorPtrPair HtslibBcfFacade::iterate(const UnpackPolicy level) const
{
    HtsBcfSrPtr sr {bcf_sr_init(), HtsSrsDeleter {}};
    
    if (bcf_sr_add_reader(sr.get(), file_path_.c_str()) != 1) {
        sr.release();
        throw std::runtime_error {"failed to open file " + file_path_.string()};
    }
    
    return std::make_pair(std::make_unique<RecordIterator>(*this, std::move(sr), level),
                          std::make_unique<RecordIterator>(*this));
}

HtslibBcfFacade::RecordIteratorPtrPair
HtslibBcfFacade::iterate(const std::string& contig, const UnpackPolicy level) const
{
    HtsBcfSrPtr sr {bcf_sr_init(), HtsSrsDeleter {}};
    
    if (bcf_sr_set_regions(sr.get(), contig.c_str(), 0) != 0) {
        throw std::runtime_error {"failed load contig " + contig};
    }
    if (bcf_sr_add_reader(sr.get(), file_path_.c_str()) != 1) {
        sr.release();
        throw std::runtime_error {"failed to open file " + file_path_.string()};
    }
    
    return std::make_pair(std::make_unique<RecordIterator>(*this, std::move(sr), level),
                          std::make_unique<RecordIterator>(*this));
}

HtslibBcfFacade::RecordIteratorPtrPair
HtslibBcfFacade::iterate(const GenomicRegion& region, const UnpackPolicy level) const
{
    HtsBcfSrPtr sr {bcf_sr_init(), HtsSrsDeleter {}};
    const auto region_str = to_string(region);
    
    if (bcf_sr_set_regions(sr.get(), region_str.c_str(), 0) != 0) {
        throw std::runtime_error {"failed load region " + region_str};
    }
    if (bcf_sr_add_reader(sr.get(), file_path_.c_str()) != 1) {
        sr.release();
        throw std::runtime_error {"failed to open file " + file_path_.string()};
    }
    
    return std::make_pair(std::make_unique<RecordIterator>(*this, std::move(sr), level),
                          std::make_unique<RecordIterator>(*this));
}

HtslibBcfFacade::RecordContainer
HtslibBcfFacade::fetch_records(const UnpackPolicy level) const
{
    const auto n_records = count_records();
    HtsBcfSrPtr sr {bcf_sr_init(), HtsSrsDeleter {}};
    
    if (bcf_sr_add_reader(sr.get(), file_path_.c_str()) != 1) {
        sr.release();
        throw std::runtime_error {"failed to open file " + file_path_.string()};
    }
    
    return fetch_records(sr.get(), level, n_records);
}

HtslibBcfFacade::RecordContainer
HtslibBcfFacade::fetch_records(const std::string& contig, const UnpackPolicy level) const
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
    
    return fetch_records(sr.get(), level, n_records);
}

HtslibBcfFacade::RecordContainer
HtslibBcfFacade::fetch_records(const GenomicRegion& region, const UnpackPolicy level) const
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
    
    return fetch_records(sr.get(), level, n_records);
}

auto hts_tag_type(const std::string& tag)
{
    using namespace vcfspec::header::meta::tag;
    const static std::unordered_map<std::string, int> types {
        {filter, BCF_HL_FLT},
        {info,   BCF_HL_INFO},
        {format, BCF_HL_FMT},
        {contig, BCF_HL_CTG}
    };
    return (types.count(tag) == 1) ? types.at(tag) : BCF_HL_STR;
}

void HtslibBcfFacade::write(const VcfHeader& header)
{
    if (file_ == nullptr) {
        throw std::runtime_error {"HtslibBcfFacade: trying to write header to closed file"};
    }
    
    auto hdr = bcf_hdr_init("w");
    
    bcf_hdr_set_version(hdr, header.file_format().c_str());
    
    for (auto& p : header.basic_fields()) {
        auto hrec = (bcf_hrec_t*) std::malloc(sizeof(bcf_hrec_t));
        hrec->type  = BCF_HL_GEN;
        hrec->key   = convert(p.first);
        hrec->value = convert(p.second);
        hrec->nkeys = 0;
        hrec->keys  = nullptr;
        hrec->vals  = nullptr;
        bcf_hdr_add_hrec(hdr, hrec);
    }
    
    for (auto& tag : header.tags()) {
        const auto type = hts_tag_type(tag);
        
        for (auto fields : header.structured_fields(tag)) {
            auto hrec = (bcf_hrec_t*) std::malloc(sizeof(bcf_hrec_t));
            
            hrec->type  = type;
            hrec->key   = convert(tag);
            hrec->nkeys = static_cast<int>(fields.size());
            hrec->keys  = (char**) std::malloc(sizeof(char*) * fields.size());
            hrec->vals  = (char**) std::malloc(sizeof(char*) * fields.size());
            
            unsigned i {0};
            // Make sure the reserved fields are written in the required order
            for (const auto& field : vcfspec::header::meta::struc::order) {
                if (fields.count(field) == 1) {
                    hrec->keys[i] = convert(field);
                    hrec->vals[i] = convert(fields[field]);
                    fields.erase(field);
                    ++i;
                }
            }
            // The rest of the fields go in whatever order they come in the map
            for (const auto& field : fields) {
                hrec->keys[i] = convert(field.first);
                hrec->vals[i] = convert(field.second);
                ++i;
            }
            hrec->value = nullptr;
            bcf_hdr_add_hrec(hdr, hrec);
        }
    }
    
    for (const auto& sample : header.samples()) {
        bcf_hdr_add_sample(hdr, sample.c_str());
    }
    bcf_hdr_write(file_.get(), hdr);
    header_.reset(hdr);
    samples_ = extract_samples(header_.get());
}

void set_chrom(const bcf_hdr_t* header, bcf1_t* record, const std::string& chrom);
void set_pos(bcf1_t* record, GenomicRegion::Position pos);
void set_id(bcf1_t* record, const std::string& id);
void set_alleles(const bcf_hdr_t* header, bcf1_t* record, const VcfRecord::NucleotideSequence& ref,
                 const std::vector<VcfRecord::NucleotideSequence>& alts);
void set_qual(bcf1_t* record, VcfRecord::QualityType qual);
void set_filter(const bcf_hdr_t* header, bcf1_t* record, const std::vector<std::string>& filters);
void set_info(const bcf_hdr_t* header, bcf1_t* dest, const VcfRecord& source);
void set_samples(const bcf_hdr_t* header, bcf1_t* dest, const VcfRecord& source,
                 const std::vector<std::string>& samples);

void HtslibBcfFacade::write(const VcfRecord& record)
{
    if (file_ == nullptr) {
        throw std::runtime_error {"HtslibBcfFacade: trying to write record to closed file"};
    }
    if (header_ == nullptr) {
        throw std::runtime_error {"HtslibBcfFacade: trying to write record without a header"};
    }
    
    const auto& contig = record.chrom();
    
    if (bcf_hdr_get_hrec(header_.get(), BCF_HL_CTG, "ID", contig.c_str(), nullptr) == nullptr) {
        throw std::runtime_error {"HtslibBcfFacade: required contig header line missing for contig \"" + contig + "\""};
    }
    
    auto hts_record = bcf_init();
    set_chrom(header_.get(), hts_record, contig);
    set_pos(hts_record, record.pos() - 1);
    set_id(hts_record, record.id());
    set_alleles(header_.get(), hts_record, record.ref(), record.alt());
    if (record.qual()) {
        set_qual(hts_record, *record.qual());
    }
    set_filter(header_.get(), hts_record, record.filter());
    set_info(header_.get(), hts_record, record);
    if (record.num_samples() > 0) {
        set_samples(header_.get(), hts_record, record, samples_);
    }
    bcf_write(file_.get(), header_.get(), hts_record);
    bcf_destroy(hts_record);
}

// HtslibBcfFacade::RecordIterator

HtslibBcfFacade::RecordIterator::RecordIterator(const HtslibBcfFacade& facade)
: facade_ {facade}
, hts_iterator_ {nullptr}
, level_ {}
, record_ {nullptr}
{}

HtslibBcfFacade::RecordIterator::RecordIterator(const HtslibBcfFacade& facade,
                                                HtsBcfSrPtr hts_iterator,
                                                UnpackPolicy level)
: facade_ {facade}
, hts_iterator_ {std::move(hts_iterator)}
, level_ {level}
{
    if (bcf_sr_next_line(hts_iterator_.get())) {
        record_ = std::make_shared<VcfRecord>(facade_.get().fetch_record(hts_iterator_.get(), level_));
    } else {
        hts_iterator_ = nullptr;
    }
}

HtslibBcfFacade::RecordIterator::reference HtslibBcfFacade::RecordIterator::operator*() const
{
    return *record_;
}

HtslibBcfFacade::RecordIterator::pointer HtslibBcfFacade::RecordIterator::operator->() const
{
    return record_.get();
}

void HtslibBcfFacade::RecordIterator::next()
{
    if (bcf_sr_next_line(hts_iterator_.get())) {
        *record_ = facade_.get().fetch_record(hts_iterator_.get(), level_);
    } else {
        hts_iterator_ = nullptr;
    }
}

HtslibBcfFacade::RecordIterator& HtslibBcfFacade::RecordIterator::operator++()
{
    this->next();
    return *this;
}

bool operator==(const HtslibBcfFacade::RecordIterator& lhs, const HtslibBcfFacade::RecordIterator& rhs)
{
    return lhs.hts_iterator_ == rhs.hts_iterator_;
}

bool operator!=(const HtslibBcfFacade::RecordIterator& lhs, const HtslibBcfFacade::RecordIterator& rhs)
{
    return !operator==(lhs, rhs);
}

// private and non-member methods

std::unordered_map<std::string, std::string> extract_format(const bcf_hrec_t* line)
{
    std::unordered_map<std::string, std::string> result {};
    result.reserve(line->nkeys);
    for (decltype(line->nkeys) k {0}; k < line->nkeys; ++k) {
        if (std::strcmp(line->keys[k], "IDX") != 0) {
            result.emplace(line->keys[k], line->vals[k]);
        }
    }
    return result;
}

void extract_chrom(const bcf_hdr_t* header, const bcf1_t* record, VcfRecord::Builder& builder)
{
    builder.set_chrom(bcf_hdr_id2name(header, record->rid));
}

void set_chrom(const bcf_hdr_t* header, bcf1_t* record, const std::string& chrom)
{
    record->rid = bcf_hdr_name2id(header, chrom.c_str());
}

void extract_pos(const bcf1_t* record, VcfRecord::Builder& builder)
{
    builder.set_pos(record->pos + 1);
}

void set_pos(bcf1_t* record, const GenomicRegion::Position pos)
{
    record->pos = static_cast<std::uint32_t>(pos);
}

void extract_id(const bcf1_t* record, VcfRecord::Builder& builder)
{
    builder.set_id(record->d.id);
}

void set_id(bcf1_t* record, const std::string& id)
{
    record->d.id = convert(id);
}

void extract_ref(const bcf1_t* record, VcfRecord::Builder& builder)
{
    builder.set_ref(record->d.allele[0]);
}

void set_alleles(const bcf_hdr_t* header, bcf1_t* record, const VcfRecord::NucleotideSequence& ref,
                 const std::vector<VcfRecord::NucleotideSequence>& alts)
{
    const auto num_alleles = alts.size() + 1;
    std::vector<const char*> alleles(num_alleles);
    alleles.front() = ref.c_str();
    std::transform(std::begin(alts), std::end(alts), std::next(std::begin(alleles)),
                   [] (const auto& allele) { return allele.c_str(); });
    bcf_update_alleles(header, record, alleles.data(), static_cast<int>(num_alleles));
}

void extract_alt(const bcf1_t* record, VcfRecord::Builder& builder)
{
    const auto num_alleles = record->n_allele;
    std::vector<VcfRecord::NucleotideSequence> alleles {};
    alleles.reserve(num_alleles - 1); // first is the reference
    for (unsigned i {1}; i < num_alleles; ++i) {
        alleles.emplace_back(record->d.allele[i]);
    }
    builder.set_alt(std::move(alleles));
}

void extract_qual(const bcf1_t* record, VcfRecord::Builder& builder)
{
    if (!std::isnan(record->qual)) {
        builder.set_qual(record->qual);
    }
}

void set_qual(bcf1_t* record, const VcfRecord::QualityType qual)
{
    record->qual = (qual == -0) ? 0 : static_cast<float>(qual);
}

void extract_filter(const bcf_hdr_t* header, const bcf1_t* record, VcfRecord::Builder& builder)
{
    std::vector<VcfRecord::KeyType> filter {};
    filter.reserve(record->d.n_flt);
    for (decltype(record->d.n_flt) i {0}; i < record->d.n_flt; ++i) {
        filter.emplace_back(bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[i]));
    }
    builder.set_filter(std::move(filter));
}

void set_filter(const bcf_hdr_t* header, bcf1_t* record, const std::vector<std::string>& filters)
{
    for (const auto& filter : filters) {
        bcf_add_filter(header, record, bcf_hdr_id2int(header, BCF_DT_ID, filter.c_str()));
    }
}

void extract_info(const bcf_hdr_t* header, bcf1_t* record, VcfRecord::Builder& builder)
{
    int* intinfo {nullptr};
    float* floatinfo {nullptr};
    char* stringinfo {nullptr};
    int* flaginfo {nullptr}; // not actually populated
    
    builder.reserve_info(record->n_info);
    
    for (unsigned i {0}; i < record->n_info; ++i) {
        int nintinfo {0};
        int nfloatinfo {0};
        int nstringinfo {0};
        int nflaginfo {0};
        const auto key_id = record->d.info[i].key;
        
        if (key_id >= header->n[BCF_DT_ID]) {
            throw std::runtime_error {"HtslibBcfFacade: found INFO key not present in header file"};
        }
        
        const char* key {header->id[BCF_DT_ID][key_id].key};
        std::vector<std::string> values {};
        
        switch (bcf_hdr_id2type(header, BCF_HL_INFO, key_id)) {
            case BCF_HT_INT:
                if (bcf_get_info_int32(header, record, key, &intinfo, &nintinfo) > 0) {
                    values.reserve(nintinfo);
                    std::transform(intinfo, intinfo + nintinfo, std::back_inserter(values),
                                   [] (auto v) {
                                       return v != bcf_int32_missing ? std::to_string(v) : vcfMissingValue;
                                   });
                }
                break;
            case BCF_HT_REAL:
                if (bcf_get_info_float(header, record, key, &floatinfo, &nfloatinfo) > 0) {
                    values.reserve(nfloatinfo);
                    std::transform(floatinfo, floatinfo + nfloatinfo, std::back_inserter(values),
                                   [] (auto v) {
                                       return v != bcf_float_missing ? std::to_string(v) : vcfMissingValue;
                                   });
                }
                break;
            case BCF_HT_STR:
            {
                const auto nchars = bcf_get_info_string(header, record, key, &stringinfo, &nstringinfo);
                if (nchars > 0) {
                    std::string tmp(stringinfo, nchars);
                    values = utils::split(tmp, vcfspec::info::valueSeperator);
                }
                break;
            }
            case BCF_HT_FLAG:
                values.reserve(1);
                values.emplace_back((bcf_get_info_flag(header, record, key, &flaginfo, &nflaginfo) == 1) ? "1" : "0");
                break;
        }
        
        builder.set_info(key, std::move(values));
    }
    
    if (intinfo != nullptr) std::free(intinfo);
    if (floatinfo != nullptr) std::free(floatinfo);
    if (stringinfo != nullptr) std::free(stringinfo);
    if (flaginfo != nullptr) std::free(flaginfo);
}

void set_info(const bcf_hdr_t* header, bcf1_t* dest, const VcfRecord& source)
{
    for (const auto& key : source.info_keys()) {
        const auto& values    = source.info_value(key);
        const auto num_values = static_cast<int>(values.size());
        
        static constexpr std::size_t defaultBufferCapacity {100};
        
        switch (bcf_hdr_id2type(header, BCF_HL_INFO, bcf_hdr_id2int(header, BCF_DT_ID, key.c_str()))) {
            case BCF_HT_INT:
            {
                bc::small_vector<int, defaultBufferCapacity> vals(num_values);
                std::transform(std::cbegin(values), std::cend(values), std::begin(vals),
                               [] (const auto& v) {
                                   return v != vcfMissingValue ? std::stoi(v) : bcf_int32_missing;
                               });
                bcf_update_info_int32(header, dest, key.c_str(), vals.data(), num_values);
                break;
            }
            case BCF_HT_REAL:
            {
                bc::small_vector<float, defaultBufferCapacity> vals(num_values);
                std::transform(std::cbegin(values), std::cend(values), std::begin(vals),
                               [] (const auto& v) {
                                   return v != vcfMissingValue ? std::stof(v) : bcf_float_missing;
                               });
                bcf_update_info_float(header, dest, key.c_str(), vals.data(), num_values);
                break;
            }
            case BCF_HT_STR:
            {
                // Can we also use small_vector here?
                const auto vals = utils::join(values, vcfspec::info::valueSeperator);
                bcf_update_info_string(header, dest, key.c_str(), vals.c_str());
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

auto extract_format(const bcf_hdr_t* header, const bcf1_t* record)
{
    std::vector<VcfRecord::KeyType> result {};
    result.reserve(record->n_fmt);
    
    for (unsigned i {0}; i < record->n_fmt; ++i) {
        const auto key_id = record->d.fmt[i].id;
        if (key_id >= header->n[BCF_DT_ID]) {
            throw std::runtime_error {"HtslibBcfFacade: found FORMAT key not present in header file"};
        }
        result.emplace_back(header->id[BCF_DT_ID][key_id].key);
    }
    
    return result;
}

void extract_samples(const bcf_hdr_t* header, bcf1_t* record, VcfRecord::Builder& builder)
{
    auto format = extract_format(header, record);
    const auto num_samples = record->n_sample;
    builder.reserve_samples(num_samples);
    
    if (format.front() == vcfspec::format::genotype) { // the first key must be GT if present
        int ngt {}, g {};
        int* gt {nullptr};
        bcf_get_genotypes(header, record, &gt, &ngt); // mallocs gt
        const auto max_ploidy = static_cast<unsigned>(record->d.fmt->n);
        
        for (unsigned sample {0}, i {0}; sample < num_samples; ++sample, i += max_ploidy) {
            std::vector<VcfRecord::NucleotideSequence> alleles {};
            alleles.reserve(max_ploidy);
            
            for (unsigned p {0}; p < max_ploidy; ++p) {
                g = gt[i + p];
                if (g == bcf_int32_vector_end) {
                    alleles.shrink_to_fit();
                    break;
                } else if (bcf_gt_is_missing(g)) {
                    alleles.push_back(vcfMissingValue);
                } else {
                    const auto idx = bcf_gt_allele(g);
                    if (idx < record->n_allele) {
                        alleles.emplace_back(record->d.allele[idx]);
                    } else {
                        alleles.push_back(vcfMissingValue);
                    }
                }
            }
            
            using Phasing = VcfRecord::Builder::Phasing;
            builder.set_genotype(header->samples[sample], std::move(alleles),
                                 bcf_gt_is_phased(g) ? Phasing::phased : Phasing::unphased);
        }
        
        std::free(gt);
    }
    
    int nintformat {};
    int* intformat {nullptr};
    int nfloatformat {};
    float* floatformat {nullptr};
    int nstringformat {};
    char** stringformat {nullptr};
    
    for (auto it = std::next(std::cbegin(format)), end = std::cend(format); it != end; ++it) {
        const auto& key = *it;
        
        std::vector<std::vector<std::string>> values(num_samples, std::vector<std::string> {});
        
        switch (bcf_hdr_id2type(header, BCF_HL_FMT, bcf_hdr_id2int(header, BCF_DT_ID, key.c_str()))) {
            case BCF_HT_INT:
                if (bcf_get_format_int32(header, record, key.c_str(), &intformat, &nintformat) > 0) {
                    const auto num_values_per_sample = nintformat / num_samples;
                    auto ptr = intformat;
                    for (unsigned sample {0}; sample < num_samples; ++sample, ptr += num_values_per_sample) {
                        values[sample].reserve(num_values_per_sample);
                        std::transform(ptr, ptr + num_values_per_sample, std::back_inserter(values[sample]),
                                       [] (auto v) {
                                           return v != bcf_int32_missing ? std::to_string(v) : vcfMissingValue;
                                       });
                    }
                }
                break;
            case BCF_HT_REAL:
                if (bcf_get_format_float(header, record, key.c_str(), &floatformat, &nfloatformat) > 0) {
                    const auto num_values_per_sample = nfloatformat / num_samples;
                    auto ptr = floatformat;
                    for (unsigned sample {0}; sample < num_samples; ++sample, ptr += num_values_per_sample) {
                        values[sample].reserve(num_values_per_sample);
                        
                        std::transform(ptr, ptr + num_samples, std::back_inserter(values[sample]),
                                       [] (auto v) {
                                           return v != bcf_float_missing ? std::to_string(v) : vcfMissingValue;
                                       });
                    }
                }
                break;
            case BCF_HT_STR:
                // TODO: Check this usage is correct. What if more than one value per sample?
                if (bcf_get_format_string(header, record, key.c_str(), &stringformat, &nstringformat) > 0) {
                    unsigned sample {0};
                    std::for_each(stringformat, stringformat + num_samples,
                                  [&values, &sample] (const char* str) {
                                      values[sample++].emplace_back(str);
                                  });
                }
                break;
        }
        
        for (unsigned sample {0}; sample < num_samples; ++sample) {
            builder.set_format(header->samples[sample], key, std::move(values[sample]));
        }
    }
    
    builder.set_format(std::move(format));
    
    if (intformat != nullptr) std::free(intformat);
    if (floatformat != nullptr) std::free(floatformat);
    if (stringformat != nullptr) {
        // bcf_get_format_string allocates two arrays
        std::free(stringformat[0]);
        std::free(stringformat);
    }
}

template <typename T, typename Container>
auto genotype_number(const T& allele, const Container& alleles, const bool is_phased)
{
    if (allele == vcfMissingValue) {
        return (is_phased) ? bcf_gt_missing + 1 : bcf_gt_missing;
    }
    const auto it = std::find(std::cbegin(alleles), std::cend(alleles), allele);
    const auto allele_num = 2 * static_cast<decltype(bcf_gt_missing)>(std::distance(std::cbegin(alleles), it)) + 2;
    return (is_phased) ? allele_num + 1 : allele_num;
}

void set_samples(const bcf_hdr_t* header, bcf1_t* dest, const VcfRecord& source,
                 const std::vector<std::string>& samples)
{
    if (samples.empty()) return;
    const auto num_samples = static_cast<int>(source.num_samples());
    const auto& format = source.format();
    if (format.empty()) return;
    
    auto first_format = std::cbegin(format);
    if (*first_format == vcfspec::format::genotype) {
        const auto& alt_alleles = source.alt();
        bc::small_vector<VcfRecord::NucleotideSequence, 5> alleles {};
        alleles.reserve(alt_alleles.size() + 1);
        alleles.push_back(source.ref());
        alleles.insert(std::end(alleles), std::cbegin(alt_alleles), std::cend(alt_alleles));
        
        unsigned max_ploidy {};
        for (const auto& sample : samples) {
            const auto p = source.ploidy(sample);
            if (p > max_ploidy) max_ploidy = p;
        }
        
        const auto ngt = num_samples * static_cast<int>(max_ploidy);
        bc::small_vector<int, 1000> genotype(ngt);
        auto it = std::begin(genotype);
        
        for (const auto& sample : samples) {
            const bool is_phased {source.is_sample_phased(sample)};
            const auto& genotype = source.get_sample_value(sample, vcfspec::format::genotype);
            const auto ploidy = static_cast<unsigned>(genotype.size());
            
            it = std::transform(std::cbegin(genotype), std::cend(genotype), it,
                                [is_phased, &alleles] (const auto& allele) {
                                    return genotype_number(allele, alleles, is_phased);
                                });
            it = std::fill_n(it, max_ploidy - ploidy, bcf_int32_vector_end);
        }
        
        bcf_update_genotypes(header, dest, genotype.data(), ngt);
        ++first_format;
    }
    
    std::for_each(first_format, std::cend(format), [&] (const auto& key) {
        const auto num_values = num_samples * static_cast<int>(source.format_cardinality(key));

        static constexpr std::size_t defaultValueCapacity {100};
        
        switch (bcf_hdr_id2type(header, BCF_HL_FMT, bcf_hdr_id2int(header, BCF_DT_ID, key.c_str()))) {
          case BCF_HT_INT:
          {
              bc::small_vector<int, defaultValueCapacity> typed_values(num_values);
              auto it = std::begin(typed_values);
              for (const auto& sample : samples) {
                  const auto& values = source.get_sample_value(sample, key);
                  it = std::transform(std::cbegin(values), std::cend(values), it,
                                      [] (const auto& v) {
                                          return v != vcfMissingValue ? std::stoi(v) : bcf_int32_missing;
                                      });
              }
              bcf_update_format_int32(header, dest, key.c_str(), typed_values.data(), num_values);
              break;
          }
          case BCF_HT_REAL:
          {
              bc::small_vector<float, defaultValueCapacity> typed_values(num_values);
              auto it = std::begin(typed_values);
              for (const auto& sample : samples) {
                  const auto& values = source.get_sample_value(sample, key);
                  it = std::transform(std::cbegin(values), std::cend(values), it,
                                      [] (const auto& v) {
                                          return v != vcfMissingValue ? std::stof(v) : bcf_float_missing;
                                      });
              }
              bcf_update_format_float(header, dest, key.c_str(), typed_values.data(), num_values);
              break;
          }
          case BCF_HT_STR:
          {
              bc::small_vector<const char*, defaultValueCapacity> typed_values(num_values);
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
    });
}

std::size_t HtslibBcfFacade::count_records(HtsBcfSrPtr& sr) const
{
    std::size_t result {0};
    while (bcf_sr_next_line(sr.get())) ++result;
    return result;
}

VcfRecord HtslibBcfFacade::fetch_record(const bcf_srs_t* sr, UnpackPolicy level) const
{
    auto hts_record = bcf_sr_get_line(sr, 0);
    bcf_unpack(hts_record, level == UnpackPolicy::all ? BCF_UN_ALL : BCF_UN_SHR);
    VcfRecord::Builder record_builder {};
    extract_chrom(header_.get(), hts_record, record_builder);
    extract_pos(hts_record, record_builder);
    extract_id(hts_record, record_builder);
    extract_ref(hts_record, record_builder);
    extract_alt(hts_record, record_builder);
    extract_qual(hts_record, record_builder);
    extract_filter(header_.get(), hts_record, record_builder);
    extract_info(header_.get(), hts_record, record_builder);
    if (level == UnpackPolicy::all && has_samples(header_.get())) {
        extract_samples(header_.get(), hts_record, record_builder);
    }
    return record_builder.build_once();
}

HtslibBcfFacade::RecordContainer
HtslibBcfFacade::fetch_records(bcf_srs_t* sr, const UnpackPolicy level, const std::size_t num_records) const
{
    RecordContainer result {};
    result.reserve(num_records);
    while (bcf_sr_next_line(sr)) {
        result.push_back(fetch_record(sr, level));
    }
    return result;
}

} // namespace octopus
