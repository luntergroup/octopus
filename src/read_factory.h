//
//  read_factory.h
//  Octopus
//
//  Created by Daniel Cooke on 14/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_factory__
#define __Octopus__read_factory__

#include <string>
#include <vector>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <iterator>

#include "genomic_region.h"
#include "read_reader.h"
#include "aligned_read.h"
#include "read_reader_implementor.h"
#include "htslib_facade.h"

class ReadFactory
{
public:
    ReadFactory() = default;
    ReadFactory(std::vector<std::string>&& read_file_paths, unsigned Max_open_files = 20);
    
    ReadFactory(const ReadFactory&)            = delete;
    ReadFactory& operator=(const ReadFactory&) = delete;
    ReadFactory(ReadFactory&&)                 = default;
    ReadFactory& operator=(ReadFactory&&)      = default;
    
    unsigned get_num_samples() const;
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

inline unsigned ReadFactory::get_num_samples() const
{
    return num_samples_;
}

// TODO: Evaluate if this should change to a member variable.
inline std::vector<std::string> ReadFactory::get_sample_ids() const
{
    std::vector<std::string> result {};
    result.reserve(num_samples_);
    for (const auto& pair : files_containing_sample_) {
        result.emplace_back(pair.first);
    }
    return result;
}

inline ReadReader ReadFactory::make_read_reader(const std::string& read_file_path)
{
    return ReadReader {read_file_path, std::make_unique<HtslibFacade>(read_file_path)};
}

inline void ReadFactory::open_reader(const std::string &read_file_path)
{
    open_readers_.emplace(std::make_pair(read_file_path, make_read_reader(read_file_path)));
    closed_files_.erase(read_file_path);
}

inline void ReadFactory::close_reader(const std::string &read_file_path)
{
    open_readers_.erase(read_file_path);
    closed_files_.insert(read_file_path);
}

inline std::vector<std::string>
ReadFactory::get_files_containing_region(const GenomicRegion& a_region) const
{
    if (files_containing_region_.count(a_region) > 0) {
        return files_containing_region_.at(a_region);
    }
    // TODO: improve this
    std::vector<std::string> result {};
    for (const auto& pair : files_containing_region_) {
        if (overlaps(pair.first, a_region)) {
            result.insert(std::end(result), std::cbegin(pair.second), std::cend(pair.second));
        }
    }
    return result;
}

#endif /* defined(__Octopus__read_factory__) */
