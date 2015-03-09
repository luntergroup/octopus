//
//  read_manager.h
//  Octopus
//
//  Created by Daniel Cooke on 14/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_manager__
#define __Octopus__read_manager__

#include <string>
#include <vector>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <memory> // std::make_unique
#include <boost/filesystem/path.hpp>

#include "genomic_region.h"
#include "read_reader.h"
#include "read_reader_impl.h"
#include "hash_functions.h"

namespace fs = boost::filesystem;

class AlignedRead;

// TODO: make this thread-safe

class ReadManager
{
public:
    using SampleIdType = IReadReaderImpl::SampleIdType;
    
    ReadManager() = default;
    explicit ReadManager(std::vector<std::string> read_file_paths, unsigned Max_open_files = 20);
    
    ReadManager(const ReadManager&)            = delete;
    ReadManager& operator=(const ReadManager&) = delete;
    ReadManager(ReadManager&&)                 = default;
    ReadManager& operator=(ReadManager&&)      = default;
    
    unsigned get_num_samples() const noexcept;
    std::vector<SampleIdType> get_sample_ids() const;
    std::vector<AlignedRead> fetch_reads(const SampleIdType& a_sample_id, const GenomicRegion& a_region);
    
private:
    using OpenReaderMap         = std::unordered_map<fs::path, ReadReader>;
    using ClosedReaders         = std::unordered_set<fs::path>;
    using SampleIdToFilePathMap = std::unordered_map<SampleIdType, std::vector<fs::path>>;
    using RegionToFilePathMap   = std::unordered_map<GenomicRegion, std::vector<fs::path>>;
    
    const unsigned Max_open_files_;
    unsigned num_samples_;
    OpenReaderMap open_readers_;
    ClosedReaders closed_readers_;
    SampleIdToFilePathMap files_containing_sample_;
    RegionToFilePathMap files_containing_region_;
    
    void setup();
    ReadReader make_read_reader(const fs::path& read_file_path);
    void open_reader(const fs::path& read_file_path);
    void close_reader(const fs::path& read_file_path);
    std::vector<fs::path> get_files_containing_region(const GenomicRegion& a_region) const;
};

#endif /* defined(__Octopus__read_manager__) */
