// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_manager_hpp
#define read_manager_hpp

#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <initializer_list>
#include <cstddef>
#include <mutex>

#include <boost/filesystem.hpp>

#include "basics/contig_region.hpp"
#include "basics/genomic_region.hpp"
#include "containers/mappable_map.hpp"
#include "utils/hash_functions.hpp"
#include "read_reader.hpp"
#include "read_reader_impl.hpp"

namespace octopus {

class AlignedRead;

namespace io {

class ReadManager
{
public:
    using Path = boost::filesystem::path;
    
    using SampleName    = IReadReaderImpl::SampleName;
    using ReadContainer = IReadReaderImpl::ReadContainer;
    using SampleReadMap = IReadReaderImpl::SampleReadMap;
    using AlignedReadReadVisitor = IReadReaderImpl::AlignedReadReadVisitor;
    using ContigRegionVisitor    = IReadReaderImpl::ContigRegionVisitor;
    
    ReadManager() = default;
    
    ReadManager(std::vector<Path> read_file_paths, unsigned max_open_files);
    ReadManager(std::initializer_list<Path> read_file_paths);
    
    ReadManager(const ReadManager&)            = delete;
    ReadManager& operator=(const ReadManager&) = delete;
    ReadManager(ReadManager &&);
    ReadManager& operator=(ReadManager &&);
    
    ~ReadManager() = default;
    
    friend void swap(ReadManager& lhs, ReadManager& rhs) noexcept;
    
    void close() const noexcept; // close all readers
    bool good() const noexcept;
    unsigned num_files() const noexcept;
    std::vector<Path> paths() const; // Managed files
    bool all_readers_have_one_sample() const;
    unsigned num_samples() const noexcept;
    const std::vector<SampleName>& samples() const;
    unsigned drop_samples(std::vector<SampleName> samples);
    
    void iterate(const GenomicRegion& region,
                 AlignedReadReadVisitor visitor) const;
    void iterate(const SampleName& sample,
                 const GenomicRegion& region,
                 AlignedReadReadVisitor visitor) const;
    void iterate(const std::vector<SampleName>& samples,
                 const GenomicRegion& region,
                 AlignedReadReadVisitor visitor) const;
    
    void iterate(const GenomicRegion& region,
                 ContigRegionVisitor visitor) const;
    void iterate(const SampleName& sample,
                 const GenomicRegion& region,
                 ContigRegionVisitor visitor) const;
    void iterate(const std::vector<SampleName>& samples,
                 const GenomicRegion& region,
                 ContigRegionVisitor visitor) const;
    
    bool has_reads(const SampleName& sample, const GenomicRegion& region) const;
    bool has_reads(const std::vector<SampleName>& samples, const GenomicRegion& region) const;
    bool has_reads(const GenomicRegion& region) const;
    
    std::size_t count_reads(const SampleName& sample, const GenomicRegion& region) const;
    std::size_t count_reads(const std::vector<SampleName>& samples, const GenomicRegion& region) const;
    std::size_t count_reads(const GenomicRegion& region) const;
    
    GenomicRegion find_covered_subregion(const SampleName& sample, const GenomicRegion& region,
                                         std::size_t max_reads) const;
    GenomicRegion find_covered_subregion(const std::vector<SampleName>& samples, const GenomicRegion& region,
                                         std::size_t max_reads) const;
    GenomicRegion find_covered_subregion(const GenomicRegion& region, std::size_t max_reads) const;
    
    ReadContainer fetch_reads(const SampleName& sample,  const GenomicRegion& region) const;
    SampleReadMap fetch_reads(const std::vector<SampleName>& samples, const GenomicRegion& region) const;
    SampleReadMap fetch_reads(const GenomicRegion& region) const;
    ReadContainer fetch_reads(const SampleName& sample,  const std::vector<GenomicRegion>& regions) const;
    SampleReadMap fetch_reads(const std::vector<SampleName>& samples, const std::vector<GenomicRegion>& regions) const;
    SampleReadMap fetch_reads(const std::vector<GenomicRegion>& regions) const;
    
private:
    using PathHash = octopus::utils::FilepathHash;
    
    struct FileSizeCompare
    {
        bool operator()(const Path& lhs, const Path& rhs) const;
    };
    
    using OpenReaderMap           = std::map<Path, ReadReader, FileSizeCompare>;
    using ClosedReaderSet         = std::unordered_set<Path, PathHash>;
    using SampleIdToReaderPathMap = std::unordered_map<SampleName, std::vector<Path>>;
    using ContigMap               = MappableMap<GenomicRegion::ContigName, ContigRegion>;
    using ReaderRegionsMap        = std::unordered_map<Path, ContigMap, PathHash>;
    
    unsigned max_open_files_ = 200;
    unsigned num_files_;
    bool all_readers_single_sample_;
    
    mutable ClosedReaderSet closed_readers_;
    mutable OpenReaderMap open_readers_;
    
    SampleIdToReaderPathMap reader_paths_containing_sample_;
    ReaderRegionsMap possible_regions_in_readers_;
    std::vector<SampleName> samples_;
    
    mutable std::mutex mutex_;
    
    void setup_reader_samples_and_regions();
    void open_initial_files();
    
    ReadReader make_reader(const Path& reader_path) const;
    bool all_readers_are_open() const noexcept;
    bool is_open(const Path& reader_path) const noexcept;
    std::vector<Path>::iterator partition_open(std::vector<Path>& reader_paths) const;
    unsigned num_open_readers() const noexcept;
    unsigned num_reader_spaces() const noexcept;
    void open_reader(const Path& reader_path) const;
    std::vector<Path>::iterator open_readers(std::vector<Path>::iterator first,
                                             std::vector<Path>::iterator last) const;
    void open_readers(unsigned n) const;
    void close_reader(const Path& reader_path) const;
    Path choose_reader_to_close() const;
    void close_readers(unsigned n) const;
    
    template <typename Visitor>
    void iterate_helper(const std::vector<SampleName>& samples,
                        const GenomicRegion& region,
                        Visitor visitor) const;
    
    void add_possible_regions_to_reader_map(const Path& reader_path, const std::vector<GenomicRegion>& regions);
    void add_reader_to_sample_map(const Path& reader_path, const std::vector<SampleName>& samples_in_reader);
    bool could_reader_contain_region(const Path& reader_path, const GenomicRegion& region) const;
    
    std::vector<Path> get_reader_paths_containing_samples(const std::vector<SampleName>& sample) const;
    std::vector<Path> get_possible_reader_paths(const GenomicRegion& region) const;
    std::vector<Path> get_possible_reader_paths(const std::vector<SampleName>& samples,
                                                const GenomicRegion& region) const;
    std::vector<Path> get_possible_reader_paths(const std::vector<GenomicRegion>& regions) const;
    std::vector<Path> get_possible_reader_paths(const std::vector<SampleName>& samples,
                                                const std::vector<GenomicRegion>& regions) const;
};

template <typename Visitor>
void ReadManager::iterate_helper(const std::vector<SampleName>& samples,
                                 const GenomicRegion& region,
                                 Visitor visitor) const
{
    if (all_readers_are_open()) {
        for (const auto& p : open_readers_) {
            if (!p.second.iterate(samples, region, visitor)) return;
        }
    } else {
        std::lock_guard<std::mutex> lock {mutex_};
        auto reader_paths = get_reader_paths_containing_samples(samples);
        auto reader_itr = partition_open(reader_paths);
        while (!reader_paths.empty()) {
            auto itr = std::find_if(reader_itr, end(reader_paths), [&] (const auto& reader_path) {
                return !open_readers_.at(reader_path).iterate(samples, region, visitor);
            });
            if (itr != end(reader_paths)) break;
            reader_paths.erase(reader_itr, end(reader_paths));
            reader_itr = open_readers(begin(reader_paths), end(reader_paths));
        }
    }
}

} // namespace io

using io::ReadManager;

} // namespace octopus

#endif
