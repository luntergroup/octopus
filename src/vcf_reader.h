//
//  vcf_reader.h
//  Octopus
//
//  Created by Daniel Cooke on 28/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__vcf_reader__
#define __Octopus__vcf_reader__

#include <vector>
#include <boost/filesystem/path.hpp>

#include "htslib_bcf_facade.h"

class VcfHeader;
class VcfRecord;
class GenomicRegion;

namespace fs = boost::filesystem;

class VcfReader
{
public:
    VcfReader()  = delete;
    explicit VcfReader(const fs::path& file_path);
    ~VcfReader() = default;
    
    VcfReader(const VcfReader&)            = default;
    VcfReader& operator=(const VcfReader&) = default;
    VcfReader(VcfReader&&)                 = default;
    VcfReader& operator=(VcfReader&&)      = default;
    
    VcfHeader fetch_header();
    std::vector<VcfRecord> fetch_records();
    std::vector<VcfRecord> fetch_records(const GenomicRegion& region);
    
private:
    fs::path file_path_;
    HtslibBcfFacade reader_;
};

#endif /* defined(__Octopus__vcf_reader__) */
