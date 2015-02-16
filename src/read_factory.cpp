//
//  read_factory.cpp
//  Octopus
//
//  Created by Daniel Cooke on 14/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "read_factory.h"

#include <algorithm>

ReadFactory::ReadFactory(std::vector<std::string>&& read_file_paths, unsigned Max_open_files)
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

void ReadFactory::setup()
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

std::vector<AlignedRead>
ReadFactory::fetch_reads(const std::string& a_sample_id, const GenomicRegion& a_region)
{
    std::vector<AlignedRead> result {};
    auto good_files = get_files_containing_region(a_region);
    for (auto& file : good_files) {
        if (open_readers_.count(file) == 0) {
            open_reader(file);
        }
        auto reads = open_readers_.at(file).fetch_reads(a_region).at(a_sample_id);
        result.insert(std::end(result), std::make_move_iterator(std::begin(reads)),
                      std::make_move_iterator(std::end(reads)));
    }
    return result;
}