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
#include <cstddef>

class GenomicRegion;
class VcfHeader;
class VcfRecord;

class IVcfReaderImpl
{
public:
    enum class Unpack { All, AllButSamples };
    
    virtual VcfHeader fetch_header() = 0;
    virtual std::size_t num_records() = 0;
    virtual std::size_t num_records(const GenomicRegion& region) = 0;
    virtual std::vector<VcfRecord> fetch_records(Unpack level = Unpack::All) = 0; // fetches all records
    virtual std::vector<VcfRecord> fetch_records(const GenomicRegion& region, Unpack level = Unpack::All) = 0;
    
    virtual ~IVcfReaderImpl() noexcept = default;
};

#endif
