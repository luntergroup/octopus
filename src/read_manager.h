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
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <memory> // std::make_unique
#include <boost/filesystem.hpp>

#include "genomic_region.h"
#include "read_reader.h"
#include "read_reader_impl.h"
#include "hash_functions.h"

namespace fs = boost::filesystem;

auto file_size_compare = [] (const fs::path& lhs, const fs::path& rhs) {
    return fs::file_size(lhs) < fs::file_size(rhs);
};

class AlignedRead;

// TODO: make this thread-safe

class ReadManager
{
public:
    using SampleIdType  = IReadReaderImpl::SampleIdType;
    using SampleReadMap = std::unordered_map<SampleIdType, std::vector<AlignedRead>>;
    
    ReadManager() = default;
    explicit ReadManager(std::vector<std::string> read_file_paths, unsigned Max_open_files = 20);
    
    ReadManager(const ReadManager&)            = delete;
    ReadManager& operator=(const ReadManager&) = delete;
    ReadManager(ReadManager&&)                 = default;
    ReadManager& operator=(ReadManager&&)      = default;
    
    unsigned get_num_samples() const noexcept;
    std::vector<SampleIdType> get_sample_ids() const;
    std::vector<AlignedRead> fetch_reads(const SampleIdType& a_sample_id, const GenomicRegion& a_region);
    SampleReadMap fetch_reads(const std::vector<SampleIdType>& sample_ids, const GenomicRegion& a_region);
    
private:
    using OpenReaderMap           = std::map<fs::path, ReadReader, decltype(file_size_compare)>;
    using ClosedReaders           = std::unordered_set<fs::path>;
    using SampleIdToReaderPathMap = std::unordered_map<SampleIdType, std::vector<fs::path>>;
    using RegionToReaderPathMap   = std::unordered_map<fs::path, std::vector<GenomicRegion>>;
    
    const unsigned Max_open_files_;
    OpenReaderMap open_readers_;
    ClosedReaders closed_readers_;
    SampleIdToReaderPathMap reader_paths_containing_sample_;
    RegionToReaderPathMap possible_regions_in_readers_;
    
    void setup();
    void check_files_exists() const;
    void get_reader_samples_and_regions();
    void open_initial_files();
    
    ReadReader make_read_reader(const fs::path& a_reader_path);
    bool is_open(const fs::path& a_reader_path) const noexcept;
    void open_reader(const fs::path& a_reader_path);
    void open_readers(std::vector<fs::path>::iterator first, std::vector<fs::path>::iterator last);
    void close_reader(const fs::path& a_reader_path);
    fs::path choose_reader_to_close() const;
    void close_readers(unsigned n);
    bool reader_could_contain_region(const fs::path& the_reader_path, const GenomicRegion& a_region) const;
    std::vector<fs::path> get_reader_paths_containing_samples(const std::vector<SampleIdType>& sample_ids) const;
    std::vector<fs::path> get_reader_paths_possibly_containing_region(const GenomicRegion& a_region) const;
};

#endif /* defined(__Octopus__read_manager__) */
