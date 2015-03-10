//
//  read_manager.cpp
//  Octopus
//
//  Created by Daniel Cooke on 14/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "read_manager.h"

#include <iterator> // std::make_move_iterator
#include <algorithm> // std::move, std::any_of, std::min

#include "htslib_read_facade.h"
#include "aligned_read.h"

ReadManager::ReadManager(std::vector<std::string> read_file_paths, unsigned Max_open_files)
:   closed_readers_ (std::make_move_iterator(std::begin(read_file_paths)),
                   std::make_move_iterator(std::end(read_file_paths))),
    open_readers_ {},
    Max_open_files_ {Max_open_files},
    num_samples_ {0},
    reader_paths_containing_sample_ {},
    possible_regions_in_readers_ {}
{
    setup();
}

void ReadManager::setup()
{
    get_reader_info();
    open_initial_files();
}

void ReadManager::get_reader_info()
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
    auto num_files_to_open = std::min(Max_open_files_, static_cast<unsigned>(closed_readers_.size()));
    for (const auto& read_file : closed_readers_) {
        if (num_files_to_open == 0) break;
        open_reader(read_file);
        --num_files_to_open;
    }
}

unsigned ReadManager::get_num_samples() const noexcept
{
    return num_samples_;
}

// TODO: Evaluate if this should change to a member variable.
std::vector<ReadManager::SampleIdType> ReadManager::get_sample_ids() const
{
    std::vector<SampleIdType> result {};
    result.reserve(num_samples_);
    for (const auto& pair : reader_paths_containing_sample_) {
        result.emplace_back(pair.first);
    }
    return result;
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
    open_readers_.emplace(a_reader_path, make_read_reader(a_reader_path));
    closed_readers_.erase(a_reader_path);
}

void ReadManager::close_reader(const fs::path& read_file_path)
{
    open_readers_.erase(read_file_path);
    closed_readers_.insert(read_file_path);
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

std::vector<AlignedRead>
ReadManager::fetch_reads(const SampleIdType& a_sample_id, const GenomicRegion& a_region)
{
    std::vector<AlignedRead> result {};
    
    for (const auto& reader_path : reader_paths_containing_sample_[a_sample_id]) {
        if (!reader_could_contain_region(reader_path, a_region)) {
            continue;
        }
        
        if (!is_open(reader_path)) {
            open_reader(reader_path);
        }
        
        auto reads = open_readers_.at(reader_path).fetch_reads(a_region);
        result.insert(std::end(result),
                      std::make_move_iterator(std::begin(reads.at(a_sample_id))),
                      std::make_move_iterator(std::end(reads.at(a_sample_id))));
    }
    
    return result;
}