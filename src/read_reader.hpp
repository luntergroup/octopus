//
//  read_reader.hpp
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_reader__
#define __Octopus__read_reader__

#include <vector>
#include <cstddef> // size_t
#include <iterator>
#include <memory>  // std::unique_ptr
#include <boost/filesystem/path.hpp>

#include "read_reader_impl.hpp"
#include "genomic_region.hpp"
#include "aligned_read.hpp"
#include "equitable.hpp"

namespace fs = boost::filesystem;

class ReadReader : public Equitable<ReadReader>
{
public:
    using SampleIdType       = IReadReaderImpl::SampleIdType;
    using SizeType           = IReadReaderImpl::SizeType;
    using SampleIdToReadsMap = IReadReaderImpl::SampleIdToReadsMap;
    
    ReadReader() = delete;
    explicit ReadReader(const fs::path& file_path, std::unique_ptr<IReadReaderImpl> the_impl);
    
    ReadReader(const ReadReader&)            = delete;
    ReadReader& operator=(const ReadReader&) = delete;
    ReadReader(ReadReader&&)                 = default;
    ReadReader& operator=(ReadReader&&)      = default;
    
    void open();
    void close();
    
    const fs::path& get_read_file_path() const noexcept;
    std::vector<std::string> get_reference_contig_names();
    std::vector<SampleIdType> get_samples();
    std::vector<std::string> get_read_groups_in_sample(const SampleIdType& sample);
    unsigned get_num_reference_contigs();
    std::vector<GenomicRegion> get_possible_regions_in_file();
    size_t count_reads(const GenomicRegion& region);
    size_t count_reads(const SampleIdType& sample, const GenomicRegion& region);
    GenomicRegion find_head_region(const GenomicRegion& region, size_t target_coverage);
    SampleIdToReadsMap fetch_reads(const GenomicRegion& region);
    std::vector<AlignedRead> fetch_reads(const SampleIdType& sample, const GenomicRegion& region);
    
private:
    fs::path file_path_;
    std::unique_ptr<IReadReaderImpl> the_impl_;
};

inline ReadReader::ReadReader(const fs::path& file_path, std::unique_ptr<IReadReaderImpl> the_impl)
:
file_path_ {file_path},
the_impl_ {std::move(the_impl)}
{}

inline void ReadReader::open()
{
    the_impl_->open();
}

inline void ReadReader::close()
{
    the_impl_->close();
}

inline const fs::path& ReadReader::get_read_file_path() const noexcept
{
    return file_path_;
}

inline std::vector<ReadReader::SampleIdType> ReadReader::get_samples()
{
    return the_impl_->get_samples();
}

inline std::vector<std::string> ReadReader::get_read_groups_in_sample(const SampleIdType& sample)
{
    return the_impl_->get_read_groups_in_sample(sample);
}

inline unsigned ReadReader::get_num_reference_contigs()
{
    return the_impl_->get_num_reference_contigs();
}

inline size_t ReadReader::count_reads(const GenomicRegion& region)
{
    return the_impl_->count_reads(region);
}

inline size_t ReadReader::count_reads(const SampleIdType& sample, const GenomicRegion& region)
{
    return the_impl_->count_reads(sample, region);
}

inline GenomicRegion ReadReader::find_head_region(const GenomicRegion& region, size_t target_coverage)
{
    return the_impl_->find_head_region(region, target_coverage);
}

inline ReadReader::SampleIdToReadsMap ReadReader::fetch_reads(const GenomicRegion& region)
{
    return the_impl_->fetch_reads(region);
}

inline std::vector<AlignedRead> ReadReader::fetch_reads(const SampleIdType& sample, const GenomicRegion& region)
{
    return the_impl_->fetch_reads(sample, region);
}

inline std::vector<std::string> ReadReader::get_reference_contig_names()
{
    return the_impl_->get_reference_contig_names();
}

inline std::vector<GenomicRegion> ReadReader::get_possible_regions_in_file()
{
    return the_impl_->get_possible_regions_in_file();
}

inline bool operator==(const ReadReader& lhs, const ReadReader& rhs)
{
    return lhs.get_read_file_path() == rhs.get_read_file_path();
}

namespace std {
    template <> struct hash<ReadReader>
    {
        size_t operator()(const ReadReader& a_read_reader) const
        {
            return hash<string>()(a_read_reader.get_read_file_path().string());
        }
    };
}

#endif /* defined(__Octopus__read_reader__) */
