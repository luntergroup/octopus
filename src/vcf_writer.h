//
//  vcf_writer.h
//  Octopus
//
//  Created by Daniel Cooke on 29/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__vcf_writer__
#define __Octopus__vcf_writer__

#include <boost/filesystem/path.hpp>

#include "htslib_bcf_facade.h"

class VcfHeader;
class VcfRecord;
class GenomicRegion;

namespace fs = boost::filesystem;

class VcfWriter
{
public:
    VcfWriter()  = delete;
    explicit VcfWriter(const fs::path& file_path);
    explicit VcfWriter(const fs::path& file_path, const VcfHeader& header);
    ~VcfWriter() = default;
    
    VcfWriter(const VcfWriter&)            = default;
    VcfWriter& operator=(const VcfWriter&) = default;
    VcfWriter(VcfWriter&&)                 = default;
    VcfWriter& operator=(VcfWriter&&)      = default;
    
    const fs::path path() const;
    void write(const VcfHeader& header);
    void write(const VcfRecord& record);
    
private:
    fs::path file_path_;
    bool is_header_written_;
    HtslibBcfFacade writer_;
};

#endif /* defined(__Octopus__vcf_writer__) */
