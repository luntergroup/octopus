//
//  vcf_reader.cpp
//  Octopus
//
//  Created by Daniel Cooke on 28/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "vcf_reader.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <utility>

#include "vcf_header.hpp"
#include "vcf_record.hpp"
#include "htslib_bcf_facade.hpp"
#include "vcf_parser.hpp"

std::unique_ptr<IVcfReaderImpl> make_vcf_reader(const VcfReader::Path& file_path)
{
    auto file_type = file_path.extension().string();
    
    if (file_type == ".vcf") {
        auto vcf_file_size = fs::file_size(file_path);
        if (vcf_file_size > 1e9) { // 1GB
            throw std::runtime_error {"VCF file " + file_path.string() + " is too big"};
        }
        
        return std::make_unique<VcfParser>(file_path);
    } else {
        return std::make_unique<HtslibBcfFacade>(file_path, "r");
    }
}

VcfReader::VcfReader(Path file_path)
:
file_path_ {std::move(file_path)},
reader_ {make_vcf_reader(file_path)}
{}

VcfReader::VcfReader(VcfReader&& other)
{
    std::lock_guard<std::mutex> lock {other.mutex_};
    file_path_ = std::move(other.file_path_);
    reader_  = std::move(other.reader_);
}

bool VcfReader::is_open() const noexcept
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_ != nullptr;
}

void VcfReader::open(Path file_path) noexcept
{
    std::lock_guard<std::mutex> lock {mutex_};
    try {
        file_path_ = std::move(file_path);
        reader_    = std::make_unique<HtslibBcfFacade>(file_path_, "w");
    } catch (...) {
        this->close();
    }
}

void VcfReader::close() noexcept
{
    std::lock_guard<std::mutex> lock {mutex_};
    reader_.reset(nullptr);
    file_path_.clear();
}

const VcfReader::Path VcfReader::path() const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return file_path_;
}

VcfHeader VcfReader::fetch_header() const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_->fetch_header();
}

size_t VcfReader::count_records()
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_->count_records();
}

size_t VcfReader::count_records(const std::string& contig)
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_->count_records(contig);
}

size_t VcfReader::count_records(const GenomicRegion& region)
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_->count_records(region);
}

std::vector<VcfRecord> VcfReader::fetch_records(Unpack level)
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_->fetch_records((level == Unpack::All) ? HtslibBcfFacade::Unpack::All : HtslibBcfFacade::Unpack::AllButSamples);
}

std::vector<VcfRecord> VcfReader::fetch_records(const std::string& contig, Unpack level)
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_->fetch_records(contig, (level == Unpack::All) ? HtslibBcfFacade::Unpack::All : HtslibBcfFacade::Unpack::AllButSamples);
}

std::vector<VcfRecord> VcfReader::fetch_records(const GenomicRegion& region, Unpack level)
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_->fetch_records(region, (level == Unpack::All) ? HtslibBcfFacade::Unpack::All : HtslibBcfFacade::Unpack::AllButSamples);
}

// non member methods

bool operator==(const VcfReader& lhs, const VcfReader& rhs)
{
    return lhs.path() == rhs.path();
}
