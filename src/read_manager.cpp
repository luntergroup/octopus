//
//  read_manager.cpp
//  Octopus
//
//  Created by Daniel Cooke on 14/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "read_manager.h"

#include <iterator> // std::make_move_iterator
#include <algorithm> // std::move

#include "htslib_facade.h"
#include "aligned_read.h"

ReadManager::ReadManager(std::vector<std::string>&& read_file_paths, unsigned Max_open_files)
:   closed_files_ (std::make_move_iterator(std::begin(read_file_paths)),
                   std::make_move_iterator(std::end(read_file_paths))),
    open_readers_ {},
    Max_open_files_ {Max_open_files},
    num_samples_ {0},
    files_containing_sample_ {},
    files_containing_region_ {}
{
    setup();
}

void ReadManager::setup()
{
    for (const auto& read_file_path : closed_files_) {
        auto read_reader = make_read_reader(read_file_path);
        for (auto&& reference_region : read_reader.get_regions_in_file()) {
            files_containing_region_[std::move(reference_region)].emplace_back(read_file_path);
        }
        for (auto&& sample_id : read_reader.get_sample_ids()) {
            files_containing_sample_[std::move(sample_id)].emplace_back(read_file_path);
        }
    }
    auto num_files_to_open = std::min(Max_open_files_, static_cast<unsigned>(closed_files_.size()));
    for (const auto& read_file : closed_files_) {
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
std::vector<std::string> ReadManager::get_sample_ids() const
{
    std::vector<std::string> result {};
    result.reserve(num_samples_);
    for (const auto& pair : files_containing_sample_) {
        result.emplace_back(pair.first);
    }
    return result;
}

ReadReader ReadManager::make_read_reader(const std::string& read_file_path)
{
    return ReadReader {read_file_path, std::make_unique<HtslibFacade>(read_file_path)};
}

void ReadManager::open_reader(const std::string &read_file_path)
{
    open_readers_.emplace(read_file_path, make_read_reader(read_file_path));
    closed_files_.erase(read_file_path);
}

void ReadManager::close_reader(const std::string &read_file_path)
{
    open_readers_.erase(read_file_path);
    closed_files_.insert(read_file_path);
}

std::vector<std::string>
ReadManager::get_files_containing_region(const GenomicRegion& a_region) const
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

std::vector<AlignedRead>
ReadManager::fetch_reads(const std::string& a_sample_id, const GenomicRegion& a_region)
{
    std::vector<AlignedRead> result {};
    for (const auto& file : get_files_containing_region(a_region)) {
        if (open_readers_.count(file) == 0) {
            open_reader(file);
        }
        auto&& reads = open_readers_.at(file).fetch_reads(a_region);
        if (reads.count(a_sample_id) > 0) {
            auto&& reads_for_sample = reads.at(a_sample_id);
            result.insert(std::end(result), std::make_move_iterator(std::begin(reads_for_sample)),
                          std::make_move_iterator(std::end(reads_for_sample)));
        }
    }
    return result;
}