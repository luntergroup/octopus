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

#include "genomic_region.h"
#include "read_reader.h"
#include "read_reader_impl.h"

class AlignedRead;

// TODO: make this thread-safe

class ReadManager
{
public:
    ReadManager() = default;
    explicit ReadManager(std::vector<std::string>&& read_file_paths, unsigned Max_open_files = 20);
    
    ReadManager(const ReadManager&)            = delete;
    ReadManager& operator=(const ReadManager&) = delete;
    ReadManager(ReadManager&&)                 = default;
    ReadManager& operator=(ReadManager&&)      = default;
    
    unsigned get_num_samples() const noexcept;
    std::vector<std::string> get_sample_ids() const;
    std::vector<AlignedRead> fetch_reads(const std::string& a_sample_id, const GenomicRegion& a_region);
    
private:
    using SampleIdToFileNameMap = std::unordered_map<std::string, std::vector<std::string>>;
    using RegionToFileNameMap   = std::unordered_map<GenomicRegion, std::vector<std::string>>;
    
    const unsigned Max_open_files_;
    unsigned num_samples_;
    SampleIdToFileNameMap files_containing_sample_;
    std::unordered_map<std::string, ReadReader> open_readers_;
    std::unordered_set<std::string> closed_files_;
    RegionToFileNameMap files_containing_region_;
    
    void setup();
    ReadReader make_read_reader(const std::string& read_file_path);
    void open_reader(const std::string& read_file_path);
    void close_reader(const std::string& read_file_path);
    std::vector<std::string> get_files_containing_region(const GenomicRegion& a_region) const;
};

#endif /* defined(__Octopus__read_manager__) */
