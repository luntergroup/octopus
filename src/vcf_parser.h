//
//  vcf_parser.h
//  Octopus
//
//  Created by Daniel Cooke on 11/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__vcf_parser__
#define __Octopus__vcf_parser__

#include <vector>
#include <cstddef> // std::size_t
#include <fstream>
#include <boost/filesystem/path.hpp>

#include "i_vcf_reader_impl.h"
#include "vcf_header.h"

namespace fs = boost::filesystem;

class GenomicRegion;
class VcfRecord;

class VcfParser : public IVcfReaderImpl
{
public:
    VcfParser() = delete;
    explicit VcfParser(const fs::path& file_path);
    
    VcfParser(const VcfParser&)            = default;
    VcfParser& operator=(const VcfParser&) = default;
    VcfParser(VcfParser&&)                 = default;
    VcfParser& operator=(VcfParser&&)      = default;
    
    VcfHeader fetch_header() override;
    std::size_t num_records() override;
    std::size_t num_records(const GenomicRegion& region) override;
    std::vector<VcfRecord> fetch_records(Unpack level = Unpack::All) override;
    std::vector<VcfRecord> fetch_records(const GenomicRegion& region, Unpack level = Unpack::All) override;
    
private:
    fs::path file_path_;
    std::ifstream file_;
    VcfHeader header_;
    const std::vector<std::string> samples_;
    const std::streampos first_record_pos_; // must go after header_!
    
    void reset_vcf();
};

#endif /* defined(__Octopus__vcf_parser__) */
