//
//  vcf_reader.hpp
//  Octopus
//
//  Created by Daniel Cooke on 28/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__vcf_reader__
#define __Octopus__vcf_reader__

#include <vector>
#include <string>
#include <cstddef>
#include <memory>
#include <mutex>

#include <boost/filesystem.hpp>

#include "i_vcf_reader_impl.hpp"

class VcfHeader;
class VcfRecord;
class GenomicRegion;

namespace fs = boost::filesystem;

class VcfReader
{
public:
    // Extracting the sample data from a VCF/BCF file can be very expensive. If this data is not
    // required, performance can be vastly improved by simply not extracting it from file.
    enum class Unpack { All, AllButSamples };
    
    VcfReader()  = delete;
    explicit VcfReader(const fs::path& file_path);
    ~VcfReader() = default;
    
    VcfReader(const VcfReader&)            = delete;
    VcfReader& operator=(const VcfReader&) = delete;
    VcfReader(VcfReader&&);
    VcfReader& operator=(VcfReader&&)      = default;
    
    const fs::path path() const;
    VcfHeader fetch_header() const;
    size_t count_records();
    size_t count_records(const std::string& contig);
    size_t count_records(const GenomicRegion& region);
    std::vector<VcfRecord> fetch_records(Unpack level = Unpack::All); // fetches all records
    std::vector<VcfRecord> fetch_records(const std::string& contig, Unpack level = Unpack::All);
    std::vector<VcfRecord> fetch_records(const GenomicRegion& region, Unpack level = Unpack::All);
    
private:
    fs::path file_path_;
    std::unique_ptr<IVcfReaderImpl> reader_;
    
    std::mutex mutex_;
};

bool operator==(const VcfReader& lhs, const VcfReader& rhs);

namespace std {
    template <> struct hash<VcfReader>
    {
        size_t operator()(const VcfReader& reader) const
        {
            return hash<string>()(reader.path().string());
        }
    };
} // namespace std

namespace std
{
    template <> struct hash<reference_wrapper<const VcfReader>>
    {
        size_t operator()(reference_wrapper<const VcfReader> reader) const
        {
            return hash<VcfReader>()(reader);
        }
    };
} // namespace std

#endif /* defined(__Octopus__vcf_reader__) */
