// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_reader_hpp
#define read_reader_hpp

#include <vector>
#include <cstddef>
#include <memory>
#include <mutex>
#include <utility>
#include <functional>

#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>

#include "concepts/equitable.hpp"
#include "read_reader_impl.hpp"

namespace octopus {

class GenomicRegion;
class AlignedRead;

namespace io {

/*
 ReadReader is a simple RAII threadsafe wrapper around a IReadReaderImpl
 */
class ReadReader : public Equitable<ReadReader>
{
public:
    using Path = boost::filesystem::path;
    
    using SampleName    = IReadReaderImpl::SampleName;
    using ReadContainer = IReadReaderImpl::ReadContainer;
    using SampleReadMap = IReadReaderImpl::SampleReadMap;
    using PositionList  = IReadReaderImpl::PositionList;
    using AlignedReadReadVisitor = IReadReaderImpl::AlignedReadReadVisitor;
    using ContigRegionVisitor    = IReadReaderImpl::ContigRegionVisitor;
    
    ReadReader() = default;
    
    ReadReader(const Path& file_path);
    
    ReadReader(const ReadReader&)            = delete;
    ReadReader& operator=(const ReadReader&) = delete;
    ReadReader(ReadReader&&);
    ReadReader& operator=(ReadReader&&)      = delete;
    
    ~ReadReader() = default;
    
    friend void swap(ReadReader& lhs, ReadReader& rhs) noexcept;
    
    bool is_open() const noexcept;
    void open();
    void close();
    
    const Path& path() const noexcept;
    
    std::vector<SampleName> extract_samples() const;
    
    std::vector<std::string> extract_read_groups(const SampleName& sample) const;
    
    std::vector<GenomicRegion::ContigName> reference_contigs() const;
    GenomicRegion::Size reference_size(const GenomicRegion::ContigName& contig) const;
    boost::optional<std::vector<GenomicRegion::ContigName>> mapped_contigs() const;
    boost::optional<std::vector<GenomicRegion>> mapped_regions() const;
    
    bool iterate(const GenomicRegion& region,
                 AlignedReadReadVisitor visitor) const;
    bool iterate(const SampleName& sample,
                 const GenomicRegion& region,
                 AlignedReadReadVisitor visitor) const;
    bool iterate(const std::vector<SampleName>& samples,
                 const GenomicRegion& region,
                 AlignedReadReadVisitor visitor) const;
    
    bool iterate(const GenomicRegion& region,
                 ContigRegionVisitor visitor) const;
    bool iterate(const SampleName& sample,
                 const GenomicRegion& region,
                 ContigRegionVisitor visitor) const;
    bool iterate(const std::vector<SampleName>& samples,
                 const GenomicRegion& region,
                 ContigRegionVisitor visitor) const;
    
    bool has_reads(const GenomicRegion& region) const;
    bool has_reads(const SampleName& sample,
                   const GenomicRegion& region) const;
    bool has_reads(const std::vector<SampleName>& samples,
                   const GenomicRegion& region) const;
    
    std::size_t count_reads(const GenomicRegion& region) const;
    std::size_t count_reads(const SampleName& sample,
                            const GenomicRegion& region) const;
    std::size_t count_reads(const std::vector<SampleName>& samples,
                            const GenomicRegion& region) const;
    
    PositionList extract_read_positions(const GenomicRegion& region,
                                        std::size_t max_coverage) const;
    PositionList extract_read_positions(const SampleName& sample,
                                        const GenomicRegion& region,
                                        std::size_t max_coverage) const;
    PositionList extract_read_positions(const std::vector<SampleName>& samples,
                                        const GenomicRegion& region,
                                        std::size_t max_coverage) const;
    
    SampleReadMap fetch_reads(const GenomicRegion& region) const;
    ReadContainer fetch_reads(const SampleName& sample,
                              const GenomicRegion& region) const;
    SampleReadMap fetch_reads(const std::vector<SampleName>& samples,
                              const GenomicRegion& region) const;
    SampleReadMap fetch_reads(const std::vector<GenomicRegion>& regions) const;
    ReadContainer fetch_reads(const SampleName& sample,
                              const std::vector<GenomicRegion>& regions) const;
    SampleReadMap fetch_reads(const std::vector<SampleName>& samples,
                              const std::vector<GenomicRegion>& regions) const;
    
private:
    Path file_path_;
    std::unique_ptr<IReadReaderImpl> impl_;
    
    mutable std::mutex mutex_;
};

bool operator==(const ReadReader& lhs, const ReadReader& rhs);

} // namespace io
} // namespace octopus

namespace std {
    template <> struct hash<octopus::io::ReadReader>
    {
        size_t operator()(const octopus::io::ReadReader& reader) const
        {
            return hash<string>()(reader.path().string());
        }
    };
    
    template <> struct hash<reference_wrapper<const octopus::io::ReadReader>>
    {
        size_t operator()(reference_wrapper<const octopus::io::ReadReader> reader) const
        {
            return hash<octopus::io::ReadReader>()(reader);
        }
    };
} // namespace std

#endif
