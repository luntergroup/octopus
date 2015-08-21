//
//  read_reader.h
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_reader__
#define __Octopus__read_reader__

#include <cstddef> // std::size_t
#include <iterator>
#include <memory>  // std::unique_ptr
#include <boost/filesystem/path.hpp>

#include "read_reader_impl.h"
#include "genomic_region.h"
#include "aligned_read.h"
#include "equitable.h"

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
    std::vector<SampleIdType> get_sample_ids();
    std::vector<std::string> get_read_groups_in_sample(const SampleIdType& a_sample_id);
    unsigned get_num_reference_contigs();
    std::vector<GenomicRegion> get_possible_regions_in_file();
    std::size_t get_num_reads(const GenomicRegion& a_region);
    SampleIdToReadsMap fetch_reads(const GenomicRegion& a_region);
    
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

inline std::vector<ReadReader::SampleIdType> ReadReader::get_sample_ids()
{
    return the_impl_->get_sample_ids();
}

inline std::vector<std::string> ReadReader::get_read_groups_in_sample(const SampleIdType& a_sample_id)
{
    return the_impl_->get_read_groups_in_sample(a_sample_id);
}

inline unsigned ReadReader::get_num_reference_contigs()
{
    return the_impl_->get_num_reference_contigs();
}

inline std::size_t ReadReader::get_num_reads(const GenomicRegion& a_region)
{
    return the_impl_->get_num_reads(a_region);
}

inline ReadReader::SampleIdToReadsMap ReadReader::fetch_reads(const GenomicRegion& a_region)
{
    return the_impl_->fetch_reads(a_region);
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
