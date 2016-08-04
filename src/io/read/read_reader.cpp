//
//  read_reader.cpp
//  octopus
//
//  Created by Daniel Cooke on 20/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "read_reader.hpp"

#include <stdexcept>
#include <sstream>

#include "htslib_sam_facade.hpp"

namespace octopus { namespace io {

class UnknownReadFileFormat : public std::runtime_error
{
public:
    UnknownReadFileFormat(boost::filesystem::path file_path)
    :
    runtime_error {"Unknown read file format"},
    file_path_ {std::move(file_path)},
    msg_ {}
    {}
    
    virtual ~UnknownReadFileFormat() = default;
    
    virtual const char* what() const noexcept override
    {
        std::ostringstream ss {};
        ss << runtime_error::what() << ": the file " << file_path_ << " is not a BAM or CRAM file";
        msg_ = ss.str();
        return msg_.c_str();
    }
private:
    boost::filesystem::path file_path_;
    mutable std::string msg_;
};

namespace {
    auto make_reader(const boost::filesystem::path& file_path)
    {
        const auto file_type = file_path.extension().string();
        
        if (file_type != ".bam" && file_path != ".cram") {
            throw UnknownReadFileFormat {file_path};
        }
        
        return std::make_unique<HtslibSamFacade>(file_path);
    }
}

ReadReader::ReadReader(const boost::filesystem::path& file_path)
:
file_path_ {file_path},
impl_ {make_reader(file_path_)}
{}

ReadReader::ReadReader(ReadReader&& other)
{
    std::lock_guard<std::mutex> lock {other.mutex_};
    file_path_ = std::move(other.file_path_);
    impl_  = std::move(other.impl_);
}

void swap(ReadReader& lhs, ReadReader& rhs) noexcept
{
    using std::swap;
    if (&lhs == &rhs) return;
    std::lock(lhs.mutex_, rhs.mutex_);
    std::lock_guard<std::mutex> lock_lhs {lhs.mutex_, std::adopt_lock}, lock_rhs {rhs.mutex_, std::adopt_lock};
    swap(lhs.file_path_, rhs.file_path_);
    swap(lhs.impl_, rhs.impl_);
}

bool ReadReader::is_open() const noexcept
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->is_open();
}

void ReadReader::open()
{
    std::lock_guard<std::mutex> lock {mutex_};
    impl_->open();
}

void ReadReader::close()
{
    std::lock_guard<std::mutex> lock {mutex_};
    impl_->close();
}

const ReadReader::Path& ReadReader::path() const noexcept
{
    return file_path_;
}

std::vector<ReadReader::SampleName> ReadReader::extract_samples() const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->extract_samples();
}

std::vector<std::string> ReadReader::extract_read_groups_in_sample(const SampleName& sample) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->extract_read_groups_in_sample(sample);
}

unsigned ReadReader::count_reference_contigs() const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->count_reference_contigs();
}

bool ReadReader::has_reads(const GenomicRegion& region) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->has_reads(region);
}

bool ReadReader::has_reads(const SampleName& sample, const GenomicRegion& region) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->has_reads(sample, region);
}

bool ReadReader::has_reads(const std::vector<SampleName>& samples,
                           const GenomicRegion& region) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->has_reads(samples, region);
}

std::size_t ReadReader::count_reads(const GenomicRegion& region) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->count_reads(region);
}

std::size_t ReadReader::count_reads(const SampleName& sample, const GenomicRegion& region) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->count_reads(sample, region);
}

std::size_t ReadReader::count_reads(const std::vector<SampleName>& samples, const GenomicRegion& region) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->count_reads(samples, region);
}

ReadReader::CoveragePair
ReadReader::find_covered_subregion(const GenomicRegion& region, std::size_t max_coverage) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->find_covered_subregion(region, max_coverage);
}

ReadReader::CoveragePair
ReadReader::find_covered_subregion(const SampleName& sample, const GenomicRegion& region,
                                   std::size_t max_coverage) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->find_covered_subregion(sample, region, max_coverage);
}

ReadReader::CoveragePair
ReadReader::find_covered_subregion(const std::vector<SampleName>& samples,
                                   const GenomicRegion& region, std::size_t max_coverage) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->find_covered_subregion(samples, region, max_coverage);
}

ReadReader::SampleReadMap ReadReader::fetch_reads(const GenomicRegion& region) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->fetch_reads(region);
}

ReadReader::ReadContainer ReadReader::fetch_reads(const SampleName& sample, const GenomicRegion& region) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->fetch_reads(sample, region);
}

ReadReader::SampleReadMap ReadReader::fetch_reads(const std::vector<SampleName>& samples,
                                                  const GenomicRegion& region) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->fetch_reads(samples, region);
}

std::vector<std::string> ReadReader::extract_reference_contig_names() const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->extract_reference_contig_names();
}

std::vector<GenomicRegion> ReadReader::extract_possible_regions_in_file() const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->extract_possible_regions_in_file();
}

bool operator==(const ReadReader& lhs, const ReadReader& rhs)
{
    return lhs.path() == rhs.path();
}

} // namespace io
} // namespace octopus
