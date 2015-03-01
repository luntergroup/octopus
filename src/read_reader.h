//
//  read_reader.h
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_reader__
#define __Octopus__read_reader__

#include <cstddef>
#include <iterator>
#include <memory>

#include "read_reader_impl.h"
#include "genomic_region.h"
#include "aligned_read.h"
#include "equitable.h"

class ReadReader : Equitable<ReadReader>
{
public:
    using SampleIdToReadsMap = std::unordered_map<std::string, std::vector<AlignedRead>>;
    
    ReadReader() = delete;
    explicit ReadReader(const std::string& read_file_path,
                        std::unique_ptr<IReadReaderImpl> the_impl);
    
    ReadReader(const ReadReader&)            = delete;
    ReadReader& operator=(const ReadReader&) = delete;
    ReadReader(ReadReader&&)                 = default;
    ReadReader& operator=(ReadReader&&)      = default;
    
    const std::string& get_read_file_path() const noexcept;
    std::vector<std::string> get_reference_contig_names();
    std::vector<std::string> get_sample_ids();
    std::vector<std::string> get_read_groups_in_sample(const std::string& a_sample_id);
    std::vector<GenomicRegion> get_regions_in_file();
    SampleIdToReadsMap fetch_reads(const GenomicRegion& a_region);
    
private:
    std::string read_file_path_;
    std::unique_ptr<IReadReaderImpl> the_impl_;
};

inline ReadReader::ReadReader(const std::string& read_file_path,
                              std::unique_ptr<IReadReaderImpl> the_impl)
:read_file_path_ {read_file_path},
 the_impl_ {std::move(the_impl)}
{}

inline const std::string& ReadReader::get_read_file_path() const noexcept
{
    return read_file_path_;
}

inline std::vector<std::string> ReadReader::get_sample_ids()
{
    return the_impl_->get_sample_ids();
}

inline std::vector<std::string> ReadReader::get_read_groups_in_sample(const std::string& a_sample_id)
{
    return the_impl_->get_read_groups_in_sample(a_sample_id);
}

inline ReadReader::SampleIdToReadsMap ReadReader::fetch_reads(const GenomicRegion& a_region)
{
    return the_impl_->fetch_reads(a_region);
}

inline std::vector<std::string> ReadReader::get_reference_contig_names()
{
    return the_impl_->get_reference_contig_names();
}

inline std::vector<GenomicRegion> ReadReader::get_regions_in_file()
{
    return the_impl_->get_regions_in_file();
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
            return hash<std::string>()(a_read_reader.get_read_file_path());
        }
    };
}

#endif /* defined(__Octopus__read_reader__) */
