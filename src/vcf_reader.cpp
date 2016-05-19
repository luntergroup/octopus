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
#include <algorithm>

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
reader_ {make_vcf_reader(file_path_)}
{}

VcfReader::VcfReader(VcfReader&& other)
{
    std::lock_guard<std::mutex> lock {other.mutex_};
    file_path_ = std::move(other.file_path_);
    reader_  = std::move(other.reader_);
}

void swap(VcfReader& lhs, VcfReader& rhs) noexcept
{
    using std::swap;
    if (&lhs == &rhs) return;
    std::lock(lhs.mutex_, rhs.mutex_);
    std::lock_guard<std::mutex> lock_lhs {lhs.mutex_, std::adopt_lock}, lock_rhs {rhs.mutex_, std::adopt_lock};
    swap(lhs.file_path_, rhs.file_path_);
    swap(lhs.reader_, rhs.reader_);
}

bool VcfReader::is_open() const noexcept
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_ != nullptr;
}

void VcfReader::open() noexcept
{
    std::lock_guard<std::mutex> lock {mutex_};
    try {
        reader_ = make_vcf_reader(file_path_);
    } catch (...) {
        this->close();
    }
}

void VcfReader::close() noexcept
{
    std::lock_guard<std::mutex> lock {mutex_};
    reader_.reset(nullptr);
}

const VcfReader::Path& VcfReader::path() const noexcept
{
    return file_path_;
}

VcfHeader VcfReader::fetch_header() const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_->fetch_header();
}

std::size_t VcfReader::count_records()
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_->count_records();
}

std::size_t VcfReader::count_records(const std::string& contig)
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_->count_records(contig);
}

std::size_t VcfReader::count_records(const GenomicRegion& region)
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_->count_records(region);
}

std::vector<VcfRecord> VcfReader::fetch_records(const UnpackPolicy level)
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_->fetch_records(level);
}

std::vector<VcfRecord> VcfReader::fetch_records(const std::string& contig, const UnpackPolicy level)
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_->fetch_records(contig, level);
}

std::vector<VcfRecord> VcfReader::fetch_records(const GenomicRegion& region, const UnpackPolicy level)
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_->fetch_records(region, level);
}

// non member methods

bool operator==(const VcfReader& lhs, const VcfReader& rhs)
{
    return lhs.path() == rhs.path();
}
