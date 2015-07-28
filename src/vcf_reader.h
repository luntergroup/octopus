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

class VcfHeader;
class VcfRecord;
class GenomicRegion;

class VcfReader
{
public:
    VcfReader()  = delete;
    ~VcfReader() = default;
    
    VcfReader(const VcfReader&)            = default;
    VcfReader& operator=(const VcfReader&) = default;
    VcfReader(VcfReader&&)                 = default;
    VcfReader& operator=(VcfReader&&)      = default;
    
    VcfHeader fetch_header();
    std::vector<VcfRecord> fetch_records();
    std::vector<VcfRecord> fetch_records(const GenomicRegion& region);
private:
    
};

#endif /* defined(__Octopus__vcf_reader__) */
