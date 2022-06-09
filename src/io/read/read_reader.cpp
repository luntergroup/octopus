// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_reader.hpp"

#include <stdexcept>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>

#include "htslib_sam_facade.hpp"
#include "exceptions/user_error.hpp"

namespace octopus { namespace io {

class UnknownReadFileFormat : public UserError
{
    virtual std::string do_why() const override
    {
        std::ostringstream ss {};
        ss << "the file you specified "
           << file_path_ << " is not a BAM or CRAM file";
        return ss.str();
    }
    
    virtual std::string do_where() const override
    {
        return "make_reader";
    }
    
    virtual std::string do_help() const override
    {
        return "Enter a valid BAM or CRAM file";
    }
    
    boost::filesystem::path file_path_;
    
public:
    UnknownReadFileFormat(boost::filesystem::path file_path)
    : file_path_ {std::move(file_path)}
    {}
    
    virtual ~UnknownReadFileFormat() = default;
};

namespace {

auto get_extension(const boost::filesystem::path& file_path)
{
    std::string result {file_path.extension().string()};
    if (!result.empty() && result.front() == '.') {
        result.erase(0, 1);
    }
    return result;
}

template <typename Container>
bool includes(const Container& values, const typename Container::value_type& value)
{
    return std::find(std::cbegin(values), std::cend(values), value) != std::cend(values);
}

bool is_valid_read_file_type(const boost::filesystem::path& file_path)
{
    static const std::vector<std::string> validReadFileExtensions {"bam", "cram"};
    return includes(validReadFileExtensions, get_extension(file_path));
}

auto make_reader(const boost::filesystem::path& file_path)
{
    if (!is_valid_read_file_type(file_path)) {
        throw UnknownReadFileFormat {file_path};
    }
    return std::make_unique<HtslibSamFacade>(file_path);
}

} //namespace

ReadReader::ReadReader(const boost::filesystem::path& file_path)
: file_path_ {file_path}
, impl_ {make_reader(file_path_)}
{}

ReadReader::ReadReader(ReadReader&& other)
{
    std::lock_guard<std::mutex> lock {other.mutex_};
    file_path_ = std::move(other.file_path_);
    impl_  = std::move(other.impl_);
}

void swap(ReadReader& lhs, ReadReader& rhs) noexcept
{
    if (&lhs == &rhs) return;
    std::lock(lhs.mutex_, rhs.mutex_);
    std::lock_guard<std::mutex> lock_lhs {lhs.mutex_, std::adopt_lock}, lock_rhs {rhs.mutex_, std::adopt_lock};
    using std::swap;
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

std::vector<std::string> ReadReader::extract_read_groups(const SampleName& sample) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->extract_read_groups(sample);
}

std::vector<GenomicRegion::ContigName> ReadReader::reference_contigs() const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->reference_contigs();
}

GenomicRegion::Size ReadReader::reference_size(const GenomicRegion::ContigName& contig) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->reference_size(contig);
}

boost::optional<std::vector<GenomicRegion::ContigName>> ReadReader::mapped_contigs() const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->mapped_contigs();
}

boost::optional<std::vector<GenomicRegion>> ReadReader::mapped_regions() const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->mapped_regions();
}

bool ReadReader::iterate(const GenomicRegion& region,
                         AlignedReadReadVisitor visitor) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->iterate(region, visitor);
}

bool ReadReader::iterate(const SampleName& sample,
                         const GenomicRegion& region,
                         AlignedReadReadVisitor visitor) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->iterate(sample, region, visitor);
}

bool ReadReader::iterate(const std::vector<SampleName>& samples,
                         const GenomicRegion& region,
                         AlignedReadReadVisitor visitor) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->iterate(samples, region, visitor);
}

bool ReadReader::iterate(const GenomicRegion& region,
                         ContigRegionVisitor visitor) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->iterate(region, visitor);
}

bool ReadReader::iterate(const SampleName& sample,
                         const GenomicRegion& region,
                         ContigRegionVisitor visitor) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->iterate(sample, region, visitor);
}

bool ReadReader::iterate(const std::vector<SampleName>& samples,
                         const GenomicRegion& region,
                         ContigRegionVisitor visitor) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->iterate(samples, region, visitor);
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

ReadReader::PositionList
ReadReader::extract_read_positions(const GenomicRegion& region, std::size_t max_coverage) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->extract_read_positions(region, max_coverage);
}

ReadReader::PositionList
ReadReader::extract_read_positions(const SampleName& sample, const GenomicRegion& region,
                                   std::size_t max_coverage) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->extract_read_positions(sample, region, max_coverage);
}

ReadReader::PositionList
ReadReader::extract_read_positions(const std::vector<SampleName>& samples,
                                   const GenomicRegion& region, std::size_t max_coverage) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->extract_read_positions(samples, region, max_coverage);
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

ReadReader::SampleReadMap ReadReader::fetch_reads(const std::vector<GenomicRegion>& regions) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->fetch_reads(regions);
}

ReadReader::ReadContainer ReadReader::fetch_reads(const SampleName& sample, const std::vector<GenomicRegion>& regions) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->fetch_reads(sample, regions);
}

ReadReader::SampleReadMap ReadReader::fetch_reads(const std::vector<SampleName>& samples,
                                                  const std::vector<GenomicRegion>& regions) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return impl_->fetch_reads(samples, regions);
}

bool operator==(const ReadReader& lhs, const ReadReader& rhs)
{
    return lhs.path() == rhs.path();
}

} // namespace io
} // namespace octopus
