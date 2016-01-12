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
#include <memory>
#include <cstddef>

#include <boost/filesystem/path.hpp>

#include "i_vcf_reader_impl.hpp"

#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/synced_bcf_reader.h"

class GenomicRegion;
class VcfHeader;
class VcfRecord;

class HtslibBcfFacade : public IVcfReaderImpl
{
public:
    using Path = boost::filesystem::path;
    
    HtslibBcfFacade() = delete;
    explicit HtslibBcfFacade(Path file_path, const std::string& mode = "r");
    ~HtslibBcfFacade() noexcept override = default;
    
    HtslibBcfFacade(const HtslibBcfFacade&)            = delete;
    HtslibBcfFacade& operator=(const HtslibBcfFacade&) = delete;
    HtslibBcfFacade(HtslibBcfFacade&&)                 = default;
    HtslibBcfFacade& operator=(HtslibBcfFacade&&)      = default;
    
    VcfHeader fetch_header() const override;
    size_t count_records() override;
    size_t count_records(const std::string& contig) override;
    size_t count_records(const GenomicRegion& region) override;
    std::vector<VcfRecord> fetch_records(Unpack level) override;
    std::vector<VcfRecord> fetch_records(const std::string& contig, Unpack level) override;
    std::vector<VcfRecord> fetch_records(const GenomicRegion& region, Unpack level) override;
    
    void write_header(const VcfHeader& header);
    void write_record(const VcfRecord& record);
    
private:
    struct HtsFileDeleter
    {
        void operator()(htsFile* file) const { hts_close(file); }
    };
    struct HtsHeaderDeleter
    {
        void operator()(bcf_hdr_t* header) const { bcf_hdr_destroy(header); }
    };
    struct HtsSrsDeleter
    {
        void operator()(bcf_srs_t* file) const { bcf_sr_destroy(file); }
    };
    struct HtsBcf1Deleter
    {
        void operator()(bcf1_t* bcf1) const { bcf_destroy(bcf1); }
    };
    
    using HtsBcfSrPtr = std::unique_ptr<bcf_srs_t, HtsSrsDeleter>;
    using HtsBcf1Ptr  = std::unique_ptr<bcf1_t, HtsBcf1Deleter>;
    
    Path file_path_;
    
    std::unique_ptr<htsFile, HtsFileDeleter> file_;
    std::unique_ptr<bcf_hdr_t, HtsHeaderDeleter> header_;
    
    std::vector<std::string> samples_;
    
    size_t num_records(HtsBcfSrPtr& sr) const;
    std::vector<VcfRecord> fetch_records(HtsBcfSrPtr& sr, Unpack level, size_t num_records = 0);
};

#endif /* defined(__Octopus__htslib_bcf_facade__) */
