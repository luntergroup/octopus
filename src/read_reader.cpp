//
//  read_reader.cpp
//  Octopus
//
//  Created by Daniel Cooke on 20/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "read_reader.hpp"

#include "htslib_sam_facade.hpp"

ReadReader::ReadReader(const boost::filesystem::path& file_path)
:
file_path_ {file_path},
the_impl_ {std::make_unique<HtslibSamFacade>(file_path_)}
{}

ReadReader::ReadReader(ReadReader&& other)
{
    std::lock_guard<std::mutex> lock {other.mutex_};
    file_path_ = std::move(other.file_path_);
    the_impl_  = std::move(other.the_impl_);
}

void swap(ReadReader& lhs, ReadReader& rhs) noexcept
{
    using std::swap;
    if (&lhs == &rhs) return;
    std::lock(lhs.mutex_, rhs.mutex_);
    std::lock_guard<std::mutex> lock_lhs {lhs.mutex_, std::adopt_lock}, lock_rhs {rhs.mutex_, std::adopt_lock};
    swap(lhs.file_path_, rhs.file_path_);
    swap(lhs.the_impl_, rhs.the_impl_);
}

bool ReadReader::is_open() const noexcept
{
    std::lock_guard<std::mutex> lock {mutex_};
    return the_impl_->is_open();
}

void ReadReader::open()
{
    std::lock_guard<std::mutex> lock {mutex_};
    the_impl_->open();
}

void ReadReader::close()
{
    std::lock_guard<std::mutex> lock {mutex_};
    the_impl_->close();
}

const ReadReader::Path& ReadReader::path() const noexcept
{
    return file_path_;
}

std::vector<ReadReader::SampleIdType> ReadReader::extract_samples()
{
    std::lock_guard<std::mutex> lock {mutex_};
    return the_impl_->extract_samples();
}

std::vector<std::string> ReadReader::extract_read_groups_in_sample(const SampleIdType& sample)
{
    std::lock_guard<std::mutex> lock {mutex_};
    return the_impl_->extract_read_groups_in_sample(sample);
}

unsigned ReadReader::count_reference_contigs()
{
    std::lock_guard<std::mutex> lock {mutex_};
    return the_impl_->count_reference_contigs();
}

size_t ReadReader::count_reads(const GenomicRegion& region)
{
    std::lock_guard<std::mutex> lock {mutex_};
    return the_impl_->count_reads(region);
}

size_t ReadReader::count_reads(const SampleIdType& sample, const GenomicRegion& region)
{
    std::lock_guard<std::mutex> lock {mutex_};
    return the_impl_->count_reads(sample, region);
}

GenomicRegion ReadReader::find_covered_subregion(const GenomicRegion& region, size_t target_coverage)
{
    std::lock_guard<std::mutex> lock {mutex_};
    return the_impl_->find_covered_subregion(region, target_coverage);
}

ReadReader::SampleReadMap ReadReader::fetch_reads(const GenomicRegion& region)
{
    std::lock_guard<std::mutex> lock {mutex_};
    return the_impl_->fetch_reads(region);
}

ReadReader::Reads ReadReader::fetch_reads(const SampleIdType& sample, const GenomicRegion& region)
{
    std::lock_guard<std::mutex> lock {mutex_};
    return the_impl_->fetch_reads(sample, region);
}

ReadReader::SampleReadMap ReadReader::fetch_reads(const std::vector<SampleIdType>& samples,
                                                  const GenomicRegion& region)
{
    std::lock_guard<std::mutex> lock {mutex_};
    return the_impl_->fetch_reads(samples, region);
}

std::vector<std::string> ReadReader::extract_reference_contig_names()
{
    std::lock_guard<std::mutex> lock {mutex_};
    return the_impl_->extract_reference_contig_names();
}

std::vector<GenomicRegion> ReadReader::extract_possible_regions_in_file()
{
    std::lock_guard<std::mutex> lock {mutex_};
    return the_impl_->extract_possible_regions_in_file();
}

bool operator==(const ReadReader& lhs, const ReadReader& rhs)
{
    return lhs.path() == rhs.path();
}
