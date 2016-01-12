//
//  vcf_writer.cpp
//  Octopus
//
//  Created by Daniel Cooke on 29/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "vcf_writer.hpp"

#include <stdexcept>
#include <utility>

#include "vcf_header.hpp"
#include "vcf_record.hpp"

#include <iostream> // TEST

VcfWriter::VcfWriter(Path file_path)
:
file_path_ {std::move(file_path)},
writer_ {std::make_unique<HtslibBcfFacade>(file_path_, "w")}
{}

VcfWriter::VcfWriter(Path file_path, const VcfHeader& header)
:
VcfWriter {std::move(file_path)}
{
    this->write(std::move(header));
}

VcfWriter::VcfWriter(VcfWriter&& other)
{
    std::lock_guard<std::mutex> lock {other.mutex_};
    file_path_         = std::move(other.file_path_);
    is_header_written_ = other.is_header_written_;
    writer_            = std::move(other.writer_);
}

bool VcfWriter::is_open() const noexcept
{
    std::lock_guard<std::mutex> lock {mutex_};
    return writer_ != nullptr;
}

void VcfWriter::open(Path file_path) noexcept
{
    std::lock_guard<std::mutex> lock {mutex_};
    try {
        file_path_         = std::move(file_path);
        is_header_written_ = false;
        writer_            = std::make_unique<HtslibBcfFacade>(file_path_, "w");
    } catch (...) {
        this->close();
    }
}

void VcfWriter::close() noexcept
{
    std::lock_guard<std::mutex> lock {mutex_};
    writer_.reset(nullptr);
    is_header_written_ = false;
    file_path_.clear();
}

const VcfWriter::Path VcfWriter::path() const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return file_path_;
}

void VcfWriter::write(const VcfHeader& header)
{
    std::lock_guard<std::mutex> lock {mutex_};
    writer_->write_header(header);
    is_header_written_ = true;
}

void VcfWriter::write(const VcfRecord& record)
{
    std::lock_guard<std::mutex> lock {mutex_};
    if (is_header_written_) {
        writer_->write_record(record);
    } else {
        throw std::runtime_error {"VcfWriter::write: cannot write record as header has not been written"};
    }
}

// non member methods

bool operator==(const VcfWriter& lhs, const VcfWriter& rhs)
{
    return lhs.path() == rhs.path();
}
