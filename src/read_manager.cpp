//
//  read_manager.cpp
//  Octopus
//
//  Created by Daniel Cooke on 14/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "read_manager.h"

#include <iterator>  // std::make_move_iterator, std::cbegin etc
#include <algorithm> // std::any_of, std::min, std::nth_element, std::partition, std::for_each
#include <boost/filesystem/operations.hpp>

#include "htslib_read_facade.h"
#include "aligned_read.h"

ReadManager::ReadManager(std::vector<std::string> read_file_paths, unsigned Max_open_files)
:
closed_readers_ {std::make_move_iterator(std::begin(read_file_paths)),
                 std::make_move_iterator(std::end(read_file_paths))},
open_readers_ {file_size_compare},
Max_open_files_ {Max_open_files},
reader_paths_containing_sample_ {},
possible_regions_in_readers_ {}
{
    setup();
}

unsigned ReadManager::get_num_samples() const noexcept
{
    return static_cast<unsigned>(reader_paths_containing_sample_.size());
}

// TODO: Evaluate if this should change to a member variable.
std::vector<ReadManager::SampleIdType> ReadManager::get_sample_ids() const
{
    std::vector<SampleIdType> result {};
    result.reserve(get_num_samples());
    
    for (const auto& pair : reader_paths_containing_sample_) {
        result.emplace_back(pair.first);
    }
    
    return result;
}

std::vector<AlignedRead>
ReadManager::fetch_reads(const SampleIdType& a_sample_id, const GenomicRegion& a_region)
{
    auto& reader_paths = reader_paths_containing_sample_.at(a_sample_id);
    
    auto last_region_containing_reader = std::partition(std::begin(reader_paths), std::end(reader_paths),
                   [this, &a_region] (auto& reader_path) {
                       return reader_could_contain_region(reader_path, a_region);
                   });
    
    // So we don't have to re-open already open readers
    auto last_open_reader = std::partition(std::begin(reader_paths), last_region_containing_reader,
                                           [this] (auto& reader_path) { return is_open(reader_path); });
    
    std::vector<AlignedRead> result {};
    
    std::for_each(std::begin(reader_paths), last_open_reader,
      [this, &a_sample_id, &a_region, &result] (const auto& reader_path) {
          auto reads = std::move(open_readers_.at(reader_path).fetch_reads(a_region).at(a_sample_id));
          result.insert(std::end(result),
                        std::make_move_iterator(std::begin(reads)),
                        std::make_move_iterator(std::end(reads)));
      });
    
    open_readers(last_open_reader, last_region_containing_reader);
    
    std::for_each(last_open_reader, last_region_containing_reader,
      [this, &a_sample_id, &a_region, &result] (const auto& reader_path) {
          auto reads = std::move(open_readers_.at(reader_path).fetch_reads(a_region).at(a_sample_id));
          result.insert(std::end(result),
                        std::make_move_iterator(std::begin(reads)),
                        std::make_move_iterator(std::end(reads)));
      });
    
    return result;
}

ReadManager::SampleReadMap
ReadManager::fetch_reads(const std::vector<SampleIdType>& sample_ids,
                         const GenomicRegion& a_region)
{
    auto reader_paths_containing_samples = get_reader_paths_containing_samples(sample_ids);
    
    auto last_region_containing_reader = std::partition(std::begin(reader_paths_containing_samples),
                                                        std::end(reader_paths_containing_samples),
                                [this, &a_region] (auto& reader_path) {
                                    return reader_could_contain_region(reader_path, a_region);
                                });
    
    auto last_open_reader = std::partition(std::begin(reader_paths_containing_samples),
                                           last_region_containing_reader,
                                           [this] (auto& reader_path) { return is_open(reader_path); });
    
    SampleReadMap result {};
    
    std::for_each(std::begin(reader_paths_containing_samples), last_open_reader,
      [this, &a_region, &result] (const auto& reader_path) {
          auto reads = open_readers_.at(reader_path).fetch_reads(a_region);
          for (auto& sample_reads : reads) {
              result[sample_reads.first].insert(std::end(result[sample_reads.first]),
                                                std::make_move_iterator(std::begin(sample_reads.second)),
                                                std::make_move_iterator(std::end(sample_reads.second)));
          }
      });
    
    open_readers(last_open_reader, last_region_containing_reader);
    
    std::for_each(last_open_reader, last_region_containing_reader,
      [this, &a_region, &result] (const auto& reader_path) {
          auto reads = open_readers_.at(reader_path).fetch_reads(a_region);
          for (auto& sample_reads : reads) {
              result[sample_reads.first].insert(std::end(result[sample_reads.first]),
                                                std::make_move_iterator(std::begin(sample_reads.second)),
                                                std::make_move_iterator(std::end(sample_reads.second)));
          }
      });
    
    return result;
}

void ReadManager::setup()
{
    check_files_exists();
    get_reader_samples_and_regions();
    open_initial_files();
}

void ReadManager::check_files_exists() const
{
    for (const auto& reader_path : closed_readers_) {
        if (!fs::exists(reader_path)) {
            throw std::runtime_error {"Cannot find read file " + reader_path.string()};
        }
    }
}

void ReadManager::get_reader_samples_and_regions()
{
    for (const auto& read_file_path : closed_readers_) {
        auto read_reader = make_read_reader(read_file_path);
        possible_regions_in_readers_.emplace(read_file_path, read_reader.get_possible_regions_in_file());
        for (auto sample_id : read_reader.get_sample_ids()) {
            reader_paths_containing_sample_[sample_id].emplace_back(read_file_path);
        }
    }
}

void ReadManager::open_initial_files()
{
    std::vector<fs::path> reader_paths {std::cbegin(closed_readers_), std::cend(closed_readers_)};
    auto num_files_to_open = std::min(Max_open_files_, static_cast<unsigned>(closed_readers_.size()));
    std::nth_element(reader_paths.begin(), reader_paths.begin() + num_files_to_open, reader_paths.end(),
                     file_size_compare);
    open_readers(reader_paths.begin(), reader_paths.begin() + num_files_to_open);
}

ReadReader ReadManager::make_read_reader(const fs::path& a_reader_path)
{
    return ReadReader {a_reader_path, std::make_unique<HtslibReadFacade>(a_reader_path)};
}

bool ReadManager::is_open(const fs::path& a_reader_path) const noexcept
{
    return open_readers_.count(a_reader_path) > 0;
}

void ReadManager::open_reader(const fs::path& a_reader_path)
{
    if (open_readers_.size() == Max_open_files_) {
        close_reader(choose_reader_to_close());
    }
    open_readers_.emplace(a_reader_path, make_read_reader(a_reader_path));
    closed_readers_.erase(a_reader_path);
}

void ReadManager::open_readers(std::vector<fs::path>::iterator first, std::vector<fs::path>::iterator last)
{
    unsigned num_open_reader_spaces = Max_open_files_ - static_cast<unsigned>(open_readers_.size());
    unsigned num_requested_spaces = static_cast<unsigned>(std::distance(first, last));
    unsigned num_readers_to_close = (num_requested_spaces <= num_open_reader_spaces) ?
                0 : num_requested_spaces - num_open_reader_spaces;
    close_readers(num_readers_to_close);
    std::for_each(first, last, [this] (const auto& reader_path) {
        open_readers_.emplace(reader_path, make_read_reader(reader_path));
        closed_readers_.erase(reader_path);
    });
}

void ReadManager::close_reader(const fs::path& read_file_path)
{
    open_readers_.erase(read_file_path);
    closed_readers_.insert(read_file_path);
}

fs::path ReadManager::choose_reader_to_close() const
{
    return open_readers_.begin()->first;
}

void ReadManager::close_readers(unsigned n)
{
    for (; n > 0; --n) {
        close_reader(choose_reader_to_close());
    }
}

bool ReadManager::reader_could_contain_region(const fs::path& the_reader_path,
                                              const GenomicRegion& a_region) const
{
    const auto& possible_regions = possible_regions_in_readers_.at(the_reader_path);
    return std::any_of(std::cbegin(possible_regions), std::cend(possible_regions),
                       [&a_region] (const auto& possible_region) {
                           return overlaps(a_region, possible_region);
                       });
}

std::vector<fs::path>
ReadManager::get_reader_paths_containing_samples(const std::vector<SampleIdType>& sample_ids) const
{
    std::vector<fs::path> all_readers_containg_samples {};
    
    for (const auto& sample_id : sample_ids) {
        const auto& sample_reader_paths = reader_paths_containing_sample_.at(sample_id);
        for (const auto& reader_path : sample_reader_paths) {
            if (std::find(all_readers_containg_samples.cbegin(), all_readers_containg_samples.cend(),
                          reader_path) == all_readers_containg_samples.cend()) {
                all_readers_containg_samples.emplace_back(reader_path);
            }
        }
    }
    
    return all_readers_containg_samples;
}

std::vector<fs::path>
ReadManager::get_reader_paths_possibly_containing_region(const GenomicRegion& a_region) const
{
    std::vector<fs::path> result {};
    
    for (const auto& pair : possible_regions_in_readers_) {
        for (const auto& region_possibly_in_reader : pair.second) {
            if (overlaps(region_possibly_in_reader, a_region)) {
                result.emplace_back(pair.first);
            }
        }
    }
    
    return result;
}
