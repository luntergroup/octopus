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

#include <boost/filesystem/path.hpp>

#include "i_vcf_reader_impl.hpp"
#include "vcf_header.hpp"

namespace fs = boost::filesystem;

class GenomicRegion;
class VcfRecord;

class VcfParser : public IVcfReaderImpl
{
public:
    VcfParser() = delete;
    explicit VcfParser(const fs::path& file_path);
    
    VcfParser(const VcfParser&)            = delete;
    VcfParser& operator=(const VcfParser&) = delete;
    VcfParser(VcfParser&&)                 = default;
    VcfParser& operator=(VcfParser&&)      = default;
    
    bool is_header_written() const noexcept override;
    VcfHeader fetch_header() const override;
    size_t count_records() override;
    size_t count_records(const std::string& contig) override;
    size_t count_records(const GenomicRegion& region) override;
    std::vector<VcfRecord> fetch_records(Unpack level) override;
    std::vector<VcfRecord> fetch_records(const std::string& contig, Unpack level) override;
    std::vector<VcfRecord> fetch_records(const GenomicRegion& region, Unpack level) override;
    
private:
    fs::path file_path_;
    mutable std::ifstream file_;
    VcfHeader header_;
    const std::vector<std::string> samples_;
    const std::streampos first_record_pos_; // must go after header_!
    
    void reset_vcf();
};

#endif /* defined(__Octopus__vcf_parser__) */
