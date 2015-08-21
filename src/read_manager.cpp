//
//  read_manager.cpp
//  Octopus
//
//  Created by Daniel Cooke on 14/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "read_manager.h"

#include <iterator>  // std::make_move_iterator, std::cbegin etc
#include <algorithm> // std::copy_if, std::min, std::nth_element, std::partition, std::for_each
#include <utility>   // std::move
#include <boost/filesystem/operations.hpp>

#include "htslib_sam_facade.h"
#include "aligned_read.h"
#include "mappable_algorithms.h"

#include <iostream> // TEST

ReadManager::ReadManager(std::vector<fs::path> read_file_paths, unsigned Max_open_files)
:
closed_readers_ {std::make_move_iterator(std::begin(read_file_paths)),
                 std::make_move_iterator(std::end(read_file_paths))},
open_readers_ {detail::file_size_compare},
Max_open_files_ {Max_open_files},
reader_paths_containing_sample_ {},
possible_regions_in_readers_ {}
{
    setup();
}

ReadManager::ReadManager(std::initializer_list<fs::path> read_file_paths)
:
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

std::vector<ReadManager::SampleIdType> ReadManager::get_sample_ids() const
{
    std::vector<SampleIdType> result {};
    result.reserve(num_samples());
    
    for (const auto& pair : reader_paths_containing_sample_) {
        result.emplace_back(pair.first);
    }
    
    return result;
}

std::vector<AlignedRead> ReadManager::fetch_reads(const SampleIdType& sample_id, const GenomicRegion& region)
{
    auto& reader_paths = reader_paths_containing_sample_.at(sample_id);
    
    auto last_region_containing_reader = std::partition(std::begin(reader_paths), std::end(reader_paths),
                   [this, &region] (auto& reader_path) {
                       return could_reader_contain_region(reader_path, region);
                   });
    
    // So we don't have to re-open already open readers
    auto last_open_reader = std::partition(std::begin(reader_paths), last_region_containing_reader,
                                           [this] (auto& reader_path) { return is_open(reader_path); });
    
    std::vector<AlignedRead> result {};
    
    std::for_each(std::begin(reader_paths), last_open_reader,
      [this, &sample_id, &region, &result] (const auto& reader_path) {
          auto reads = std::move(open_readers_.at(reader_path).fetch_reads(region)[sample_id]);
          result.insert(std::end(result), std::make_move_iterator(std::begin(reads)),
                        std::make_move_iterator(std::end(reads)));
      });
    
    open_readers(last_open_reader, last_region_containing_reader);
    
    std::for_each(last_open_reader, last_region_containing_reader,
      [this, &sample_id, &region, &result] (const auto& reader_path) {
          auto reads = std::move(open_readers_.at(reader_path).fetch_reads(region)[sample_id]);
          result.insert(std::end(result), std::make_move_iterator(std::begin(reads)),
                        std::make_move_iterator(std::end(reads)));
      });
    
    return result;
}

ReadManager::SampleReadMap ReadManager::fetch_reads(const std::vector<SampleIdType>& sample_ids,
                                                    const GenomicRegion& region)
{
    auto reader_paths_containing_samples = get_reader_paths_containing_samples(sample_ids);
    
    auto last_region_containing_reader = std::partition(std::begin(reader_paths_containing_samples),
                                                        std::end(reader_paths_containing_samples),
                                                        [this, &region] (auto& reader_path) {
                                                            return could_reader_contain_region(reader_path, region);
                                                        });
    
    auto last_open_reader = std::partition(std::begin(reader_paths_containing_samples),
                                           last_region_containing_reader,
                                           [this] (auto& reader_path) { return is_open(reader_path); });
    
    SampleReadMap result {};
    
    std::for_each(std::begin(reader_paths_containing_samples), last_open_reader,
      [this, &region, &result] (const auto& reader_path) {
          auto reads = open_readers_.at(reader_path).fetch_reads(region);
          for (auto& sample_reads : reads) {
              result[sample_reads.first].insert(std::end(result[sample_reads.first]),
                                                std::make_move_iterator(std::begin(sample_reads.second)),
                                                std::make_move_iterator(std::end(sample_reads.second)));
          }
      });
    
    open_readers(last_open_reader, last_region_containing_reader);
    
    std::for_each(last_open_reader, last_region_containing_reader,
      [this, &region, &result] (const auto& reader_path) {
          auto reads = open_readers_.at(reader_path).fetch_reads(region);
          for (auto& sample_reads : reads) {
              result[sample_reads.first].insert(std::end(result[sample_reads.first]),
                                                std::make_move_iterator(std::begin(sample_reads.second)),
                                                std::make_move_iterator(std::end(sample_reads.second)));
          }
      });
    
    return result;
}

ReadManager::SampleReadMap ReadManager::fetch_reads(const GenomicRegion& region)
{
    return fetch_reads(get_sample_ids(), region);
}

// Private methods

void ReadManager::setup()
{
    auto bad_paths = get_bad_paths();
    
    if (bad_paths.empty()) {
        setup_reader_samples_and_regions();
        open_initial_files();
    } else {
        throw std::runtime_error {"Could not open read files"};
    }
}

std::vector<fs::path> ReadManager::get_bad_paths() const
{
    std::vector<fs::path> result {};
    
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
        add_reader_to_sample_map(reader_path, reader.get_sample_ids());
    }
}

void ReadManager::open_initial_files()
{
    std::vector<fs::path> reader_paths {std::cbegin(closed_readers_), std::cend(closed_readers_)};
    auto num_files_to_open = std::min(Max_open_files_, static_cast<unsigned>(closed_readers_.size()));
    
    std::nth_element(reader_paths.begin(), reader_paths.begin() + num_files_to_open, reader_paths.end(),
                     detail::file_size_compare);
    
    open_readers(reader_paths.begin(), reader_paths.begin() + num_files_to_open);
}

ReadReader ReadManager::make_reader(const fs::path& reader_path)
{
    return ReadReader {reader_path, std::make_unique<HtslibSamFacade>(reader_path)};
}

bool ReadManager::is_open(const fs::path& reader_path) const noexcept
{
    return open_readers_.count(reader_path) > 0;
}

void ReadManager::open_reader(const fs::path& reader_path)
{
    if (open_readers_.size() == Max_open_files_) {
        close_reader(choose_reader_to_close());
    }
    open_readers_.emplace(reader_path, make_reader(reader_path));
    closed_readers_.erase(reader_path);
}

void ReadManager::open_readers(std::vector<fs::path>::iterator first, std::vector<fs::path>::iterator last)
{
    unsigned num_open_reader_spaces = Max_open_files_ - static_cast<unsigned>(open_readers_.size());
    unsigned num_requested_spaces = static_cast<unsigned>(std::distance(first, last));
    unsigned num_readers_to_close = (num_requested_spaces <= num_open_reader_spaces) ?
                0 : num_requested_spaces - num_open_reader_spaces;
    
    close_readers(num_readers_to_close);
    
    std::for_each(first, last, [this] (const auto& reader_path) {
        open_readers_.emplace(reader_path, make_reader(reader_path));
        closed_readers_.erase(reader_path);
    });
}

void ReadManager::close_reader(const fs::path& read_file_path)
{
    // TODO: we can make use of IReadReaderImpl::close to avoid calling deconstructor on the file.
    open_readers_.erase(read_file_path);
    closed_readers_.insert(read_file_path);
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
    for (auto sample_id : the_samples_in_reader) {
        reader_paths_containing_sample_[sample_id].emplace_back(the_reader_path);
    }
}

std::vector<fs::path> ReadManager::get_reader_paths_containing_samples(const std::vector<SampleIdType>& sample_ids) const
{
    std::unordered_set<fs::path> unique_reader_paths {};
    
    for (const auto& sample_id : sample_ids) {
        const auto& sample_reader_paths = reader_paths_containing_sample_.at(sample_id);
        for (const auto& reader_path : sample_reader_paths) {
            unique_reader_paths.emplace(reader_path);
        }
    }
    
    return std::vector<fs::path>(unique_reader_paths.begin(), unique_reader_paths.end());
}
