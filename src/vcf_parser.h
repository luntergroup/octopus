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

#include "vcf_header.h"

namespace fs = boost::filesystem;

class GenomicRegion;
class VcfRecord;

class VcfParser
{
public:
    enum class Unpack { All, AllButSamples };
    
    VcfParser() = delete;
    explicit VcfParser(const fs::path& file_path);
    ~VcfParser() = default;
    
    VcfParser(const VcfParser&)            = default;
    VcfParser& operator=(const VcfParser&) = default;
    VcfParser(VcfParser&&)                 = default;
    VcfParser& operator=(VcfParser&&)      = default;
    
    VcfHeader fetch_header();
    std::size_t num_records() const;
    std::size_t num_records(const GenomicRegion& region) const;
    std::vector<VcfRecord> fetch_records(Unpack level = Unpack::All); // fetches all records
    std::vector<VcfRecord> fetch_records(const GenomicRegion& region, Unpack level = Unpack::All);
    
private:
    fs::path file_path_;
    std::ifstream file_;
    
    VcfHeader header_;
    
    void parse();
};

#endif /* defined(__Octopus__vcf_parser__) */
