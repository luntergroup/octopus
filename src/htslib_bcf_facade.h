//
//  htslib_bcf_facade.h
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__htslib_bcf_facade__
#define __Octopus__htslib_bcf_facade__

#include <string>
#include <set>
#include <memory> // std::unique_ptr
#include <boost/filesystem/path.hpp>

#include "variant_file_reader_impl.h"
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/synced_bcf_reader.h"

class GenomicRegion;
class Variant;
class VcfRecord;

namespace fs = boost::filesystem;

auto htslib_file_deleter       = [] (htsFile* file) { hts_close(file); };
auto htslib_bcf_header_deleter = [] (bcf_hdr_t* header) { bcf_hdr_destroy(header); };
auto htslib_bcf_srs_deleter    = [] (bcf_srs_t* file) { bcf_sr_destroy(file); };
auto htslib_bcf1_deleter       = [] (bcf1_t* bcf1) { bcf_destroy(bcf1); };

class HtslibBcfFacade : public IVariantFileReaderImpl
{
public:
    HtslibBcfFacade() = delete;
    HtslibBcfFacade(const fs::path& file_path);
    ~HtslibBcfFacade() = default;
    
    HtslibBcfFacade(const HtslibBcfFacade&)            = default;
    HtslibBcfFacade& operator=(const HtslibBcfFacade&) = default;
    HtslibBcfFacade(HtslibBcfFacade&&)                 = default;
    HtslibBcfFacade& operator=(HtslibBcfFacade&&)      = default;
    
    std::vector<Variant> fetch_variants(const GenomicRegion& region) override;
    std::vector<VcfRecord> fetch_records(const GenomicRegion& region);
    
private:
    using HtsBcfSrPtr = std::unique_ptr<bcf_srs_t, decltype(htslib_bcf_srs_deleter)>;
    using HtsBcf1Ptr  = std::unique_ptr<bcf1_t, decltype(htslib_bcf1_deleter)>;
    
    fs::path file_path_;
    std::unique_ptr<htsFile, decltype(htslib_file_deleter)> file_;
    std::unique_ptr<bcf_hdr_t, decltype(htslib_bcf_header_deleter)> header_;
};

#endif /* defined(__Octopus__htslib_bcf_facade__) */
