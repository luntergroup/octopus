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
#include <iterator>

#include <boost/filesystem/path.hpp>

#include "vcf_reader_impl.hpp"
#include "vcf_record.hpp"

#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/synced_bcf_reader.h"

class GenomicRegion;
class VcfHeader;

class HtslibBcfFacade : public IVcfReaderImpl
{
public:
    using Path = boost::filesystem::path;
    
    using IVcfReaderImpl::UnpackPolicy;
    
    using IVcfReaderImpl::RecordContainer;
    
    class RecordIterator;
    
    using IVcfReaderImpl::RecordIteratorPtrPair;
    
    HtslibBcfFacade() = delete;
    
    HtslibBcfFacade(Path file_path, const std::string& mode = "r");
    
    ~HtslibBcfFacade() noexcept override = default;
    
    HtslibBcfFacade(const HtslibBcfFacade&)            = delete;
    HtslibBcfFacade& operator=(const HtslibBcfFacade&) = delete;
    HtslibBcfFacade(HtslibBcfFacade&&)                 = default;
    HtslibBcfFacade& operator=(HtslibBcfFacade&&)      = default;
    
    bool is_header_written() const noexcept override;
    
    VcfHeader fetch_header() const override;
    
    std::size_t count_records() const override;
    std::size_t count_records(const std::string& contig) const override;
    std::size_t count_records(const GenomicRegion& region) const override;
    
    RecordIteratorPtrPair iterate(UnpackPolicy level) const override;
    //RecordIteratorPtrPair iterate(const std::string& contig, UnpackPolicy level) override;
    //RecordIteratorPtrPair iterate(const GenomicRegion& region, UnpackPolicy level) override;
    
    RecordContainer fetch_records(UnpackPolicy level) const override;
    RecordContainer fetch_records(const std::string& contig, UnpackPolicy level) const override;
    RecordContainer fetch_records(const GenomicRegion& region, UnpackPolicy level) const override;
    
    void write(const VcfHeader& header);
    void write(const VcfRecord& record);
    
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
    
    std::size_t count_records(HtsBcfSrPtr& sr) const;
    
    VcfRecord fetch_record(const bcf_srs_t* sr, UnpackPolicy level) const;
    
    RecordContainer fetch_records(bcf_srs_t*, UnpackPolicy level, size_t num_records) const;
    
    friend RecordIterator;
    
public:
    class RecordIterator : public IVcfReaderImpl::RecordIterator
    {
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type        = VcfRecord;
        using difference_type   = std::ptrdiff_t;
        using pointer           = const VcfRecord*;
        using reference         = const VcfRecord&;
        
        RecordIterator(const HtslibBcfFacade& facade);
        RecordIterator(const HtslibBcfFacade& facade, HtsBcfSrPtr hts_iterator, UnpackPolicy level);
        
        ~RecordIterator() noexcept = default;
        
        RecordIterator(const RecordIterator&)            = default;
        RecordIterator& operator=(const RecordIterator&) = default;
        RecordIterator(RecordIterator&&)                 = default;
        RecordIterator& operator=(RecordIterator&&)      = default;
        
        reference operator*() const override;
        pointer operator->() const override;
        
        void next() override;
        RecordIterator& operator++();
        
    private:
        using HtsBcfSrSharedPtr = std::shared_ptr<bcf_srs_t>;
        
        std::reference_wrapper<const HtslibBcfFacade> facade_;
        
        HtsBcfSrSharedPtr hts_iterator_;
        UnpackPolicy level_;
        
        std::shared_ptr<VcfRecord> record_;
        
        friend bool operator==(const RecordIterator& lhs, const RecordIterator& rhs);
    };
};

bool operator!=(const HtslibBcfFacade::RecordIterator& lhs, const HtslibBcfFacade::RecordIterator& rhs);

#endif /* defined(__Octopus__htslib_bcf_facade__) */
