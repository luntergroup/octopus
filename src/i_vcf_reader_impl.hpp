//
//  i_vcf_reader.h
//  Octopus
//
//  Created by Daniel Cooke on 14/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_i_vcf_reader_impl_h
#define Octopus_i_vcf_reader_impl_h

#include <vector>
#include <string>
#include <cstddef>

class GenomicRegion;
class VcfHeader;
class VcfRecord;

class IVcfReaderImpl
{
public:
    enum class UnpackPolicy { All, Sites };
    
    virtual bool is_header_written() const noexcept = 0;
    
    virtual VcfHeader fetch_header() const = 0;
    
    virtual size_t count_records() = 0;
    virtual size_t count_records(const std::string& contig) = 0;
    virtual size_t count_records(const GenomicRegion& region) = 0;
    
    virtual std::vector<VcfRecord> fetch_records(UnpackPolicy level = UnpackPolicy::All) = 0; // fetches all records
    virtual std::vector<VcfRecord> fetch_records(const std::string& contig, UnpackPolicy level = UnpackPolicy::All) = 0;
    virtual std::vector<VcfRecord> fetch_records(const GenomicRegion& region, UnpackPolicy level = UnpackPolicy::All) = 0;
    
    virtual ~IVcfReaderImpl() noexcept = default;
};

#endif
