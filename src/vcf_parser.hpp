//
//  vcf_parser.hpp
//  Octopus
//
//  Created by Daniel Cooke on 11/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__vcf_parser__
#define __Octopus__vcf_parser__

#include <vector>
#include <string>
#include <cstddef>
#include <fstream>
#include <iterator>

#include <boost/filesystem/path.hpp>

#include "vcf_reader_impl.hpp"
#include "vcf_header.hpp"
#include "vcf_record.hpp"

namespace fs = boost::filesystem;

class GenomicRegion;

class VcfParser : public IVcfReaderImpl
{
public:
    using IVcfReaderImpl::RecordContainer;
    
    class RecordIterator;
    
    using IVcfReaderImpl::RecordIteratorPtrPair;
    
    VcfParser() = delete;
    
    VcfParser(const fs::path& file_path);
    
    VcfParser(const VcfParser&)            = delete;
    VcfParser& operator=(const VcfParser&) = delete;
    VcfParser(VcfParser&&)                 = default;
    VcfParser& operator=(VcfParser&&)      = default;
    
    bool is_header_written() const noexcept override;
    
    VcfHeader fetch_header() const override;
    
    std::size_t count_records() const override;
    std::size_t count_records(const std::string& contig) const override;
    std::size_t count_records(const GenomicRegion& region) const override;
    
    RecordIteratorPtrPair iterate(UnpackPolicy level) const override;
    
    RecordContainer fetch_records(UnpackPolicy level) const override;
    RecordContainer fetch_records(const std::string& contig, UnpackPolicy level) const override;
    RecordContainer fetch_records(const GenomicRegion& region, UnpackPolicy level) const override;
    
private:
    fs::path file_path_;
    
    mutable std::ifstream file_;
    
    VcfHeader header_;
    
    const std::vector<std::string> samples_;
    const std::streampos first_record_pos_; // must go after header_!
    
    void reset_vcf() const; // logically
};

class VcfParser::RecordIterator : public IVcfReaderImpl::RecordIterator
{
public:
    using iterator_category = std::input_iterator_tag;
    using value_type        = VcfRecord;
    using difference_type   = std::ptrdiff_t;
    using pointer           = const VcfRecord*;
    using reference         = const VcfRecord&;
    
    RecordIterator() = default;
    
    ~RecordIterator() noexcept = default;
    
    RecordIterator(const RecordIterator&)            = default;
    RecordIterator& operator=(const RecordIterator&) = default;
    RecordIterator(RecordIterator&&)                 = default;
    RecordIterator& operator=(RecordIterator&&)      = default;
    
    reference operator*() const override { return record_; }
    pointer operator->() const override { return &record_; }
    
    void next() override {}
    //RecordIterator& operator++();
    
    friend bool operator==(const RecordIterator& lhs, const RecordIterator& rhs)
    {
        return lhs.record_ == rhs.record_;
    }
    
private:
    VcfRecord record_;
};

#endif /* defined(__Octopus__vcf_parser__) */
