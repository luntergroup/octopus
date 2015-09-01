//
//  read_manager.h
//  Octopus
//
//  Created by Daniel Cooke on 14/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_manager__
#define __Octopus__read_manager__

#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <initializer_list>
#include <boost/filesystem.hpp>

#include "contig_region.h"
#include "genomic_region.h"
#include "read_reader.h"
#include "read_reader_impl.h"
#include "hash_functions.h"

namespace fs = boost::filesystem;

namespace detail
{
    auto file_size_compare = [] (const fs::path& lhs, const fs::path& rhs) {
        return fs::file_size(lhs) < fs::file_size(rhs);
    };
}

class AlignedRead;

// TODO: make this thread-safe

class ReadManager
{
public:
    using SampleIdType  = IReadReaderImpl::SampleIdType;
    using SampleReadMap = std::unordered_map<SampleIdType, std::vector<AlignedRead>>;
    
    ReadManager() = default;
    explicit ReadManager(std::vector<fs::path> read_file_paths, unsigned max_open_files);
    explicit ReadManager(std::initializer_list<fs::path> read_file_paths);
    ~ReadManager() = default;
    
    ReadManager(const ReadManager&)            = delete;
    ReadManager& operator=(const ReadManager&) = delete;
    ReadManager(ReadManager&&)                 = default;
    ReadManager& operator=(ReadManager&&)      = default;
    
    unsigned num_samples() const noexcept;
    std::vector<SampleIdType> get_sample_ids() const;
    std::vector<AlignedRead> fetch_reads(const SampleIdType& sample_id, const GenomicRegion& region);
    SampleReadMap fetch_reads(const std::vector<SampleIdType>& sample_ids, const GenomicRegion& region);
    SampleReadMap fetch_reads(const GenomicRegion& region);
    
private:
    using OpenReaderMap           = std::map<fs::path, ReadReader, decltype(detail::file_size_compare)>;
    using ClosedReaders           = std::unordered_set<fs::path>;
    using SampleIdToReaderPathMap = std::unordered_map<SampleIdType, std::vector<fs::path>>;
    using ReaderRegionsMap        = std::unordered_map<fs::path, std::unordered_map<GenomicRegion::StringType,
                                                                                std::vector<ContigRegion>>>;
    
    const unsigned max_open_files_ = 200;
    const unsigned num_files_;
    OpenReaderMap open_readers_;
    ClosedReaders closed_readers_;
    SampleIdToReaderPathMap reader_paths_containing_sample_;
    ReaderRegionsMap possible_regions_in_readers_;
    
    void setup();
    std::vector<fs::path> get_bad_paths() const;
    void setup_reader_samples_and_regions();
    void open_initial_files();
    
    ReadReader make_reader(const fs::path& reader_path);
    bool is_open(const fs::path& reader_path) const noexcept;
    void open_reader(const fs::path& reader_path);
    void open_readers(std::vector<fs::path>::iterator first, std::vector<fs::path>::iterator last);
    void close_reader(const fs::path& reader_path);
    fs::path choose_reader_to_close() const;
    void close_readers(unsigned n);
    
    void add_possible_regions_to_reader_map(const fs::path& the_reader_path, const std::vector<GenomicRegion>& the_regions);
    void add_reader_to_sample_map(const fs::path& the_reader_path, const std::vector<SampleIdType>& the_samples_in_reader);
    bool could_reader_contain_region(const fs::path& the_reader_path, const GenomicRegion& region) const;
    std::vector<fs::path> get_reader_paths_containing_samples(const std::vector<SampleIdType>& sample_ids) const;
    std::vector<fs::path> get_reader_paths_possibly_containing_region(const GenomicRegion& region) const;
};

#endif /* defined(__Octopus__read_manager__) */
