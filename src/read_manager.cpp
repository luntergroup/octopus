//
//  read_manager.cpp
//  Octopus
//
//  Created by Daniel Cooke on 14/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "read_manager.hpp"

#include <memory> // std::make_unique
#include <iterator>  // std::make_move_iterator, std::cbegin, std::cend, std::distance
#include <algorithm> // std::sort, std::copy_if, std::min, std::nth_element, std::partition, std::for_each
                     // std::remove_if
#include <utility>   // std::move
#include <boost/filesystem/operations.hpp>

#include "htslib_sam_facade.hpp"
#include "aligned_read.hpp"
#include "mappable_algorithms.hpp"

#include <iostream> // TEST

ReadManager::ReadManager(std::vector<fs::path> read_file_paths, unsigned max_open_files)
:
max_open_files_ {max_open_files},
num_files_ {static_cast<unsigned>(read_file_paths.size())},
closed_readers_ {std::make_move_iterator(std::begin(read_file_paths)),
                 std::make_move_iterator(std::end(read_file_paths))},
open_readers_ {detail::file_size_compare},
reader_paths_containing_sample_ {},
possible_regions_in_readers_ {}
{
    setup();
}

ReadManager::ReadManager(std::initializer_list<fs::path> read_file_paths)
:
num_files_ {static_cast<unsigned>(read_file_paths.size())},
closed_readers_ {std::begin(read_file_paths), std::end(read_file_paths)},
open_readers_ {detail::file_size_compare},
reader_paths_containing_sample_ {},
possible_regions_in_readers_ {}
{
    setup();
}

unsigned ReadManager::num_samples() const noexcept
{
    return static_cast<unsigned>(reader_paths_containing_sample_.size());
}

std::vector<ReadManager::SampleIdType> ReadManager::get_samples() const
{
    std::vector<SampleIdType> result {};
    result.reserve(num_samples());
    
    for (const auto& pair : reader_paths_containing_sample_) {
        result.emplace_back(pair.first);
    }
    
    std::sort(std::begin(result), std::end(result)); // just for consistency
    
    return result;
}

size_t ReadManager::count_reads(const SampleIdType& sample, const GenomicRegion& region)
{
    using std::begin; using std::end; using std::for_each;
    
    auto reader_paths = get_possible_reader_paths({sample}, region);
    
    auto it = partition_open(reader_paths);
    
    size_t result {};
    
    while (!reader_paths.empty()) {
        for_each(it, end(reader_paths),
                 [this, &sample, &region, &result] (const auto& reader_path) {
                     result += open_readers_.at(reader_path).count_reads(sample, region);
                 });
        
        reader_paths.erase(it, end(reader_paths));
        it = open_readers(begin(reader_paths), end(reader_paths));
    }
    
    return result;
}

size_t ReadManager::count_reads(const std::vector<SampleIdType>& samples, const GenomicRegion& region)
{
    using std::begin; using std::end; using std::for_each;
    
    auto reader_paths = get_possible_reader_paths(samples, region);
    
    auto it = partition_open(reader_paths);
    
    size_t result {};
    
    while (!reader_paths.empty()) {
        std:for_each(it, end(reader_paths),
                     [this, &samples, &region, &result] (const auto& reader_path) {
                         for (const auto& sample : samples) {
                             result += open_readers_.at(reader_path).count_reads(sample, region);
                         }
                     });
        
        reader_paths.erase(it, end(reader_paths));
        it = open_readers(begin(reader_paths), end(reader_paths));
    }
    
    return result;
}

size_t ReadManager::count_reads(const GenomicRegion& region)
{
    return count_reads(get_samples(), region);
}

GenomicRegion ReadManager::find_head_region(const std::vector<SampleIdType>& samples,
                                            const GenomicRegion& region, size_t target_coverage)
{
    using std::begin; using std::end; using std::for_each;
    
    auto reader_paths = get_possible_reader_paths(samples, region);
    
    auto it = partition_open(reader_paths);
    
    std::vector<GenomicRegion> result {};
    
    while (!reader_paths.empty()) {
        for_each(it, end(reader_paths),
                 [this, &samples, &region, &result, target_coverage] (const auto& reader_path) {
                     result.push_back(open_readers_.at(reader_path).find_head_region(region, target_coverage));
                 });
        
        reader_paths.erase(it, end(reader_paths));
        it = open_readers(begin(reader_paths), end(reader_paths));
    }
    
    return *leftmost_mappable(std::cbegin(result), std::cend(result));
}

GenomicRegion ReadManager::find_head_region(const GenomicRegion& region, size_t target_coverage)
{
    return find_head_region(get_samples(), region, target_coverage);
}

ReadManager::Reads ReadManager::fetch_reads(const SampleIdType& sample, const GenomicRegion& region)
{
    using std::begin; using std::end; using std::move; using std::make_move_iterator; using std::for_each;
    
    auto reader_paths = get_possible_reader_paths({sample}, region);
    
    auto it = partition_open(reader_paths);
    
    Reads result {};
    
    while (!reader_paths.empty()) {
        for_each(it, end(reader_paths),
                 [this, &sample, &region, &result] (const auto& reader_path) {
                     auto reads = move(open_readers_.at(reader_path).fetch_reads(sample, region));
                     result.insert(make_move_iterator(begin(reads)), make_move_iterator(end(reads)));
                 });
        
        reader_paths.erase(it, end(reader_paths));
        it = open_readers(begin(reader_paths), end(reader_paths));
    }
    
    return result;
}

ReadManager::SampleReadMap ReadManager::fetch_reads(const std::vector<SampleIdType>& samples,
                                                    const GenomicRegion& region)
{
    using std::begin; using std::end; using std::make_move_iterator; using std::for_each;
    
    auto reader_paths = get_possible_reader_paths(samples, region);
    
    auto it = partition_open(reader_paths);
    
    SampleReadMap result {};
    result.reserve(samples.size());
    
    while (!reader_paths.empty()) {
        for_each(it, end(reader_paths),
                 [this, &region, &result] (const auto& reader_path) {
                     auto reads = open_readers_.at(reader_path).fetch_reads(region);
                     for (auto& sample_reads : reads) {
                         result[sample_reads.first].insert(make_move_iterator(begin(sample_reads.second)),
                                                           make_move_iterator(end(sample_reads.second)));
                         sample_reads.second.clear();
                     }
                 });
        
        reader_paths.erase(it, end(reader_paths));
        it = open_readers(begin(reader_paths), end(reader_paths));
    }
    
    for (const auto& sample : samples) {
        if (result.count(sample) == 0) {
            result.emplace(sample, SampleReadMap::mapped_type {});
        }
    }
    
    return result;
}

ReadManager::SampleReadMap ReadManager::fetch_reads(const GenomicRegion& region)
{
    return fetch_reads(get_samples(), region);
}

// Private methods

void ReadManager::setup()
{
    auto bad_paths = get_bad_paths();
    
    if (bad_paths.empty()) {
        setup_reader_samples_and_regions();
        open_initial_files();
    } else {
        std::string error {"bad read files: \n"};
        
        for (const auto& path : bad_paths) {
            error += "\t* " + path.string() + ": does not exist\n";
        }
        
        error.erase(--error.end()); // removes last \n
        
        throw std::runtime_error {error};
    }
}

std::vector<fs::path> ReadManager::get_bad_paths() const
{
    std::vector<fs::path> result {};
    result.reserve(num_files_);
    
    std::copy_if(std::cbegin(closed_readers_), std::cend(closed_readers_), std::back_inserter(result),
                 [] (const auto& path) {
                     return !fs::exists(path);
                 });
    
    return result;
}

void ReadManager::setup_reader_samples_and_regions()
{
    for (const auto& reader_path : closed_readers_) {
        auto reader = make_reader(reader_path);
        add_possible_regions_to_reader_map(reader_path, reader.get_possible_regions_in_file());
        add_reader_to_sample_map(reader_path, reader.get_samples());
    }
}

void ReadManager::open_initial_files()
{
    using std::begin; using std::end; using std::cbegin; using std::cend;
    
    std::vector<fs::path> reader_paths {cbegin(closed_readers_), cend(closed_readers_)};
    
    auto num_files_to_open = std::min(max_open_files_, static_cast<unsigned>(closed_readers_.size()));
    
    std::nth_element(begin(reader_paths), begin(reader_paths) + num_files_to_open, end(reader_paths),
                     detail::file_size_compare);
    
    open_readers(begin(reader_paths), begin(reader_paths) + num_files_to_open);
}

ReadReader ReadManager::make_reader(const fs::path& reader_path)
{
    return ReadReader {reader_path, std::make_unique<HtslibSamFacade>(reader_path)};
}

bool ReadManager::is_open(const fs::path& reader_path) const noexcept
{
    return open_readers_.count(reader_path) == 1;
}

std::vector<fs::path>::iterator ReadManager::partition_open(std::vector<fs::path>& reader_paths) const
{
    return std::partition(std::begin(reader_paths), std::end(reader_paths),
                          [this] (auto& reader_path) { return !is_open(reader_path); });
}

unsigned ReadManager::num_open_readers() const noexcept
{
    return static_cast<unsigned>(open_readers_.size());
}

unsigned ReadManager::num_reader_spaces() const noexcept
{
    return max_open_files_ - num_open_readers();
}

void ReadManager::open_reader(const fs::path& reader_path)
{
    if (num_open_readers() == max_open_files_) { // do we need this?
        close_reader(choose_reader_to_close());
    }
    open_readers_.emplace(reader_path, make_reader(reader_path));
    closed_readers_.erase(reader_path);
}

std::vector<fs::path>::iterator
ReadManager::open_readers(std::vector<fs::path>::iterator first, std::vector<fs::path>::iterator last)
{
    if (first == last) return first;
    
    auto num_available_spaces = num_reader_spaces();
    auto num_requested_spaces = static_cast<unsigned>(std::distance(first, last));
    
    if (num_requested_spaces <= num_available_spaces) {
        std::for_each(first, last, [this] (const auto& path) { open_reader(path); });
        return first;
    }
    
    auto num_readers_to_close = std::min(num_open_readers(), num_requested_spaces - num_available_spaces);
    
    close_readers(num_readers_to_close);
    
    num_available_spaces += num_readers_to_close;
    
    // partition range so opened readers come last
    auto first_open = std::next(first, num_requested_spaces - num_available_spaces);
    
    std::for_each(first_open, last, [this] (const auto& path) { open_reader(path); });
    
    return first_open;
}

void ReadManager::close_reader(const fs::path& reader_path)
{
    // TODO: we can make use of IReadReaderImpl::close to avoid calling deconstructor on the file.
    open_readers_.erase(reader_path);
    closed_readers_.insert(reader_path);
}

fs::path ReadManager::choose_reader_to_close() const
{
    return open_readers_.begin()->first; // i.e. smallest file size
}

void ReadManager::close_readers(unsigned n)
{
    for (; n > 0; --n) {
        close_reader(choose_reader_to_close());
    }
}

void ReadManager::add_possible_regions_to_reader_map(const fs::path& the_reader_path,
                                                     const std::vector<GenomicRegion>& the_regions)
{
    for (const auto& region : the_regions) {
        possible_regions_in_readers_[the_reader_path][region.get_contig_name()].emplace_back(region.get_contig_region());
    }
    
    for (auto& contig_regions : possible_regions_in_readers_[the_reader_path]) {
        std::sort(std::begin(contig_regions.second), std::end(contig_regions.second));
    }
}

bool ReadManager::could_reader_contain_region(const fs::path& the_reader_path, const GenomicRegion& region) const
{
    if (possible_regions_in_readers_.count(the_reader_path) == 0) return false;
    if (possible_regions_in_readers_.at(the_reader_path).count(region.get_contig_name()) == 0) return false;
    
    const auto& contig_regions = possible_regions_in_readers_.at(the_reader_path).at(region.get_contig_name());
    
    return has_overlapped(std::cbegin(contig_regions), std::cend(contig_regions), region.get_contig_region());
}

std::vector<fs::path> ReadManager::get_reader_paths_possibly_containing_region(const GenomicRegion& region) const
{
    std::vector<fs::path> result {};
    result.reserve(num_files_);
    
    for (const auto& reader_path : closed_readers_) {
        if (could_reader_contain_region(reader_path, region)) {
            result.emplace_back(reader_path);
        }
    }
    
    for (const auto& reader : open_readers_) {
        if (could_reader_contain_region(reader.first, region)) {
            result.emplace_back(reader.first);
        }
    }
    
    return result;
}

void ReadManager::add_reader_to_sample_map(const fs::path& the_reader_path,
                                           const std::vector<SampleIdType>& the_samples_in_reader)
{
    for (auto sample : the_samples_in_reader) {
        reader_paths_containing_sample_[sample].emplace_back(the_reader_path);
    }
}

std::vector<fs::path> ReadManager::get_reader_paths_containing_samples(const std::vector<SampleIdType>& samples) const
{
    std::unordered_set<fs::path> unique_reader_paths {};
    unique_reader_paths.reserve(num_files_);
    
    for (const auto& sample : samples) {
        for (const auto& reader_path : reader_paths_containing_sample_.at(sample)) {
            unique_reader_paths.emplace(reader_path);
        }
    }
    
    return std::vector<fs::path> {std::begin(unique_reader_paths), std::end(unique_reader_paths)};
}

std::vector<fs::path>
ReadManager::get_possible_reader_paths(const std::vector<SampleIdType>& samples, const GenomicRegion& region) const
{
    auto result = get_reader_paths_containing_samples(samples);
    
    auto it = std::remove_if(std::begin(result), std::end(result),
                             [this, &region] (const auto& path) { return !could_reader_contain_region(path, region); });
    
    result.erase(it, std::end(result));
    
    return result;
}
