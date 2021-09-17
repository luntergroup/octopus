// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "vcf_reader.hpp"

#include <fstream>
#include <iostream>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <typeinfo>

#include <boost/filesystem/operations.hpp>

#include "vcf_header.hpp"
#include "vcf_record.hpp"
#include "htslib_bcf_facade.hpp"
#include "vcf_parser.hpp"

namespace octopus {

std::unique_ptr<IVcfReaderImpl> make_vcf_reader(const VcfReader::Path& file_path)
{
    if (!boost::filesystem::exists(file_path)) {
        std::ostringstream ss {};
        ss << "VcfReader: the path ";
        ss << file_path;
        ss << " does not exist";
        throw std::runtime_error {ss.str()};
    }
    auto file_type = file_path.extension().string();
    if (file_type == ".vcf") {
        return std::make_unique<VcfParser>(file_path);
    } else {
        return std::make_unique<HtslibBcfFacade>(file_path, HtslibBcfFacade::Mode::read);
    }
}

// VcfReader::Iterator

template <typename T>
static decltype(auto) element_typeid(const std::unique_ptr<T>& ptr) noexcept
{
    return typeid(*ptr);
}

VcfReader::RecordIterator::RecordIterator(IVcfReaderImpl::RecordIteratorPtr itr)
: itr_ {std::move(itr)}
, type_ {element_typeid(itr_)}
{}

VcfReader::RecordIterator::RecordIterator(const RecordIterator& other)
: itr_ {other.itr_->clone()}
, type_ {other.type_}
{}

VcfReader::RecordIterator::reference VcfReader::RecordIterator::operator*() const
{
    return itr_->operator*();
}

VcfReader::RecordIterator::pointer VcfReader::RecordIterator::operator->() const
{
    return itr_->operator->();
}

VcfReader::RecordIterator& VcfReader::RecordIterator::operator++()
{
    itr_->next();
    return *this;
}

bool operator==(const VcfReader::RecordIterator& lhs, const VcfReader::RecordIterator& rhs)
{
    const static std::type_index htsType {typeid(HtslibBcfFacade::RecordIterator)};
    const static std::type_index parserType {typeid(VcfParser::RecordIterator)};
    
    if (lhs.type_ == htsType) {
        try {
            const auto& typed_lhs = dynamic_cast<HtslibBcfFacade::RecordIterator&>(*lhs.itr_);
            const auto& typed_rhs = dynamic_cast<HtslibBcfFacade::RecordIterator&>(*rhs.itr_);
            return typed_lhs == typed_rhs;
        } catch(const std::bad_cast&) {
            throw std::runtime_error {"VcfReader: trying to compare incompatible iterators"};
        }
    } else if (lhs.type_ == parserType) {
        try {
            const auto& typed_lhs = dynamic_cast<VcfParser::RecordIterator&>(*lhs.itr_);
            const auto& typed_rhs = dynamic_cast<VcfParser::RecordIterator&>(*rhs.itr_);
            return typed_lhs == typed_rhs;
        } catch (const std::bad_cast&) {
            throw std::runtime_error {"VcfReader: trying to compare incompatible iterators"};
        }
    } else {
        throw std::runtime_error {"VcfReader: trying to compare unknown iterator types"};
    }
}

VcfReader::VcfReader(Path file_path)
: file_path_ {std::move(file_path)}
, reader_ {make_vcf_reader(file_path_)}
{}

VcfReader::VcfReader(Path file_path, const ReferenceGenome& reference)
: VcfReader {std::move(file_path)}
{
    reader_->set_reference(reference);
}

VcfReader::VcfReader(VcfReader&& other)
{
    std::lock_guard<std::mutex> lock {other.mutex_};
    file_path_ = std::move(other.file_path_);
    reader_  = std::move(other.reader_);
}

VcfReader& VcfReader::operator=(VcfReader&& other)
{
    if (this != &other) {
        std::unique_lock<std::mutex> lock_lhs {mutex_, std::defer_lock}, lock_rhs {other.mutex_, std::defer_lock};
        std::lock(lock_lhs, lock_rhs);
        file_path_ = std::move(other.file_path_);
        reader_    = std::move(other.reader_);
    }
    return *this;
}

void swap(VcfReader& lhs, VcfReader& rhs) noexcept
{
    if (&lhs == &rhs) return;
    std::lock(lhs.mutex_, rhs.mutex_);
    std::lock_guard<std::mutex> lock_lhs {lhs.mutex_, std::adopt_lock}, lock_rhs {rhs.mutex_, std::adopt_lock};
    using std::swap;
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
    reader_ = make_vcf_reader(file_path_);
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

std::size_t VcfReader::count_records() const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_->count_records();
}

std::size_t VcfReader::count_records(const std::string& contig) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_->count_records(contig);
}

std::size_t VcfReader::count_records(const GenomicRegion& region) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_->count_records(region);
}

VcfReader::RecordContainer VcfReader::fetch_records(const UnpackPolicy level) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_->fetch_records(level);
}

VcfReader::RecordContainer VcfReader::fetch_records(const std::string& contig, const UnpackPolicy level) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_->fetch_records(contig, level);
}

VcfReader::RecordContainer VcfReader::fetch_records(const GenomicRegion& region, const UnpackPolicy level) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return reader_->fetch_records(region, level);
}

VcfReader::RecordIteratorPair VcfReader::iterate(const UnpackPolicy level) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    auto p = reader_->iterate(level);
    return std::make_pair(std::move(p.first), std::move(p.second));
}

VcfReader::RecordIteratorPair VcfReader::iterate(const std::string& contig, const UnpackPolicy level) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    auto p = reader_->iterate(contig, level);
    return std::make_pair(std::move(p.first), std::move(p.second));
}

VcfReader::RecordIteratorPair VcfReader::iterate(const GenomicRegion& region, const UnpackPolicy level) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    auto p = reader_->iterate(region, level);
    return std::make_pair(std::move(p.first), std::move(p.second));
}

// non member methods

bool operator==(const VcfReader& lhs, const VcfReader& rhs)
{
    return lhs.path() == rhs.path();
}

bool operator!=(const VcfReader::RecordIterator& lhs, const VcfReader::RecordIterator& rhs)
{
    return !(lhs == rhs);
}

} // namespace octopus
