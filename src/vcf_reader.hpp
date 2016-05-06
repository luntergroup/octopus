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

class VcfReader
{
public:
    using Path = boost::filesystem::path;
    
    // Extracting the sample data from a VCF/BCF file can be very expensive. If this data is not
    // required, performance can be vastly improved by simply not extracting it from file.
    enum class Unpack { All, AllButSamples };
    
    VcfReader()  = default;
    explicit VcfReader(Path file_path);
    ~VcfReader() = default;
    
    VcfReader(const VcfReader&)            = delete;
    VcfReader& operator=(const VcfReader&) = delete;
    VcfReader(VcfReader&&);
    VcfReader& operator=(VcfReader&&)      = default;
    
    friend void swap(VcfReader& lhs, VcfReader& rhs) noexcept;
    
    bool is_open() const noexcept;
    void open(Path file_path) noexcept;
    void close() noexcept;
    
    const Path path() const;
    
    VcfHeader fetch_header() const;
    
    std::size_t count_records();
    std::size_t count_records(const std::string& contig);
    std::size_t count_records(const GenomicRegion& region);
    
    std::vector<VcfRecord> fetch_records(Unpack level = Unpack::All); // fetches all records
    std::vector<VcfRecord> fetch_records(const std::string& contig, Unpack level = Unpack::All);
    std::vector<VcfRecord> fetch_records(const GenomicRegion& region, Unpack level = Unpack::All);
    
private:
    Path file_path_;
    std::unique_ptr<IVcfReaderImpl> reader_;
    
    mutable std::mutex mutex_;
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
