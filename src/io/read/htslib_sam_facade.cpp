// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "htslib_sam_facade.hpp"

#include <sstream>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <sstream>
#include <limits>
#include <cassert>
#include <array>

#include <boost/filesystem/operations.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <htslib/sam.h>

#include "basics/cigar_string.hpp"
#include "basics/genomic_region.hpp"
#include "basics/contig_region.hpp"
#include "exceptions/missing_file_error.hpp"
#include "exceptions/missing_index_error.hpp"
#include "exceptions/malformed_file_error.hpp"
#include "exceptions/unwritable_file_error.hpp"
#include "utils/string_utils.hpp"

#include <iostream>

namespace octopus { namespace io {

static const std::string readGroupTag       {"RG"};
static const std::string readGroupIdTag     {"ID"};
static const std::string sampleIdTag        {"SM"};

class MissingBAM : public MissingFileError
{
    std::string do_where() const override { return "HtslibSamFacade"; }
public:
    MissingBAM(boost::filesystem::path file) : MissingFileError {std::move(file), "bam"} {}
};


class MissingCRAM : public MissingFileError
{
    std::string do_where() const override { return "HtslibSamFacade"; }
public:
    MissingCRAM(boost::filesystem::path file) : MissingFileError {std::move(file), "cram"} {}
};

class MissingBAMIndex : public MissingIndexError
{
    std::string do_where() const override { return "HtslibSamFacade"; }
    
    std::string do_help() const override
    {
        return "ensure that a valid bam index (.bai) exists in the same directory as the given "
        "bam file. You can make one with the 'samtools index' command";
    }
public:
    MissingBAMIndex(boost::filesystem::path file) : MissingIndexError {std::move(file), "bam"} {}
};

class MissingCRAMIndex : public MissingIndexError
{
    std::string do_where() const override { return "HtslibSamFacade"; }
    
    std::string do_help() const override
    {
        return "ensure that a valid cram index (.crai) exists in the same directory as the given "
        "cram file. You can make one with the 'samtools index' command";
    }
public:
    MissingCRAMIndex(boost::filesystem::path file) : MissingIndexError {std::move(file), "cram"} {}
};

class MalformedBAM : public MalformedFileError
{
    std::string do_where() const override { return "HtslibSamFacade"; }
    
    std::string do_help() const override
    {
        return "refer to the latest SAM specification";
    }
public:
    MalformedBAM(boost::filesystem::path file) : MalformedFileError {std::move(file)} {}
};

class MalformedCRAM : public MalformedFileError
{
    std::string do_where() const override { return "HtslibSamFacade"; }
    
    std::string do_help() const override
    {
        return "refer to the latest SAM specification";
    }
public:
    MalformedCRAM(boost::filesystem::path file) : MalformedFileError {std::move(file)} {}
};

class MalformedBAMHeader : public MalformedFileError
{
    std::string do_where() const override { return "HtslibSamFacade"; }
    
    std::string do_help() const override
    {
        return "ensure the bam header contains @RG and @SM tags. Refer to the latest SAM specification for details";
    }
public:
    MalformedBAMHeader(boost::filesystem::path file) : MalformedFileError {std::move(file)} {}
};

class MalformedCRAMHeader : public MalformedFileError
{
    std::string do_where() const override { return "HtslibSamFacade"; }
    
    std::string do_help() const override
    {
        return "ensure the cram header contains @RG and @SM tags. Refer to the latest SAM specification for details";
    }
public:
    MalformedCRAMHeader(boost::filesystem::path file) : MalformedFileError {std::move(file)} {}
};

class UnwritableBAM : public UnwritableFileError
{
    std::string do_where() const override { return "HtslibSamFacade"; }
public:
    UnwritableBAM(boost::filesystem::path file) : UnwritableFileError {std::move(file), "bam"} {}
};

class InvalidBamRecord : public std::runtime_error
{
public:
    InvalidBamRecord(boost::filesystem::path file_path, std::string read_name, std::string message)
    : std::runtime_error {"Invalid BAM record"}
    , file_path_ {std::move(file_path)}
    , read_name_ {std::move(read_name)}
    , message_ {std::move(message)}
    , msg_ {}
    {}
    
    virtual ~InvalidBamRecord() = default;
    
    virtual const char* what() const noexcept override
    {
        std::ostringstream ss {};
        ss << std::runtime_error::what() << ": in file " << file_path_ << ", in read " << read_name_
                    << ": " << message_;
        msg_ = ss.str();
        return msg_.c_str();
    }
private:
    boost::filesystem::path file_path_;
    std::string read_name_, message_;
    mutable std::string msg_;
};

// public methods

namespace {

auto open_hts_file(const boost::filesystem::path& file)
{
    hts_verbose = 0; // disable hts error reporting
    return sam_open(file.c_str(), "r");
}

bool is_cram(const boost::filesystem::path& file)
{
    return file.extension().string() == ".cram";
}

} // namespace

HtslibSamFacade::HtslibSamFacade(Path file_path)
: file_path_ {std::move(file_path)}
, hts_file_ {open_hts_file(file_path_), HtsFileDeleter {}}
, hts_header_ {(hts_file_) ? sam_hdr_read(hts_file_.get()) : nullptr, HtsHeaderDeleter {}}
, hts_index_ {(hts_file_) ? sam_index_load(hts_file_.get(), file_path_.c_str()) : nullptr, HtsIndexDeleter {}}
, hts_targets_ {}
, contig_names_ {}
, sample_names_ {}
, samples_ {}
{
    namespace fs = boost::filesystem;
    if (!hts_file_) {
        if (!fs::exists(file_path_)) {
            if (is_cram(file_path_)) {
                throw MissingCRAM {file_path_};
            } else {
                throw MissingBAM {file_path_};
            }
        } else {
            if (is_cram(file_path_)) {
                throw MalformedCRAM {file_path_};
            } else {
                throw MalformedBAM {file_path_};
            }
        }
    }
    if (hts_file_ && !hts_index_) {
        if (hts_file_->is_cram) {
            throw MissingCRAMIndex {file_path};
        } else {
            throw MissingBAMIndex {file_path_};
        }
    }
    try {
        init_maps();
    } catch(...) {
        close();
        throw;
    }
    for (const auto& pair : sample_names_) {
        if (std::find(std::cbegin(samples_), std::cend(samples_), pair.second) == std::cend(samples_)) {
            samples_.emplace_back(pair.second);
        }
    }
    samples_.shrink_to_fit();
    std::sort(std::begin(samples_), std::end(samples_));
}

auto open_hts_writable_file(const boost::filesystem::path& path)
{
    std::string mode {"[w]"};
    const auto extension = path.extension();
    if (extension == ".bam") {
        mode += "b";
    } else if (extension == "cram") {
        mode += "c";
    }
    return sam_open(path.c_str(), mode.c_str());
}

HtslibSamFacade::HtslibSamFacade(Path sam_out, Path sam_template)
: HtslibSamFacade {std::move(sam_template)}
{
    file_path_ = std::move(sam_out);
    hts_file_.reset(open_hts_writable_file(file_path_));
    if (!hts_file_) {
        throw UnwritableBAM {std::move(file_path_)};
    }
    hts_index_ = nullptr;
    if (sam_hdr_write(hts_file_.get(), hts_header_.get()) < 0) {
        throw UnwritableBAM {std::move(file_path_)};
    }
}

HtslibSamFacade::~HtslibSamFacade()
{
    if (!hts_index_) {
        hts_header_.reset(nullptr);
        hts_file_.reset(nullptr);
        if (sam_index_build(file_path_.c_str(), 0) < 0) {
            return;
        }
    }
}

bool HtslibSamFacade::is_open() const noexcept
{
    return hts_file_ != nullptr && hts_header_ != nullptr && hts_index_ != nullptr;
}

void HtslibSamFacade::open()
{
    hts_file_.reset(sam_open(file_path_.string().c_str(), "r"));
    if (hts_file_) {
        hts_header_.reset(sam_hdr_read(hts_file_.get()));
        hts_index_.reset(sam_index_load(hts_file_.get(), file_path_.c_str()));
    }
}

void HtslibSamFacade::close()
{
    hts_file_.reset(nullptr);
    hts_header_.reset(nullptr);
    hts_index_.reset(nullptr);
}

GenomicRegion::Size HtslibSamFacade::reference_size(const GenomicRegion::ContigName& contig) const
{
    return hts_header_->target_len[get_htslib_target(contig)];
}

std::uint64_t HtslibSamFacade::get_num_mapped_reads(const GenomicRegion::ContigName& contig) const
{
    std::uint64_t num_mapped {}, num_unmapped {};
    hts_idx_get_stat(hts_index_.get(), get_htslib_target(contig), &num_mapped, &num_unmapped);
    return num_mapped;
}

std::vector<HtslibSamFacade::SampleName> HtslibSamFacade::extract_samples() const
{
    return samples_;
}

std::vector<std::string> HtslibSamFacade::extract_read_groups(const SampleName& sample) const
{
    std::vector<std::string> result {};
    for (const auto& pair : sample_names_) {
        if (pair.second == sample) result.emplace_back(pair.first);
    }
    result.shrink_to_fit();
    return result;
}

namespace {
    
template <typename T>
bool is_sorted(const std::vector<T>& v)
{
    return std::is_sorted(std::cbegin(v), std::cend(v));
}

template <typename T>
std::vector<T> sort(std::vector<T> v)
{
    std::sort(std::begin(v), std::end(v));
    return v;
}

template <typename T>
bool set_equal(const std::vector<T>& lhs, const std::vector<T>& rhs)
{
    if (is_sorted(lhs)) {
        return is_sorted(rhs) ? lhs == rhs : lhs == sort(rhs);
    } else if (is_sorted(rhs)) {
        return sort(lhs) == rhs;
    } else {
        return sort(lhs) == sort(rhs);
    }
}

template <typename S>
auto get_readable_samples(std::vector<S> request_samples, const std::vector<S>& file_samples)
{
    std::sort(std::begin(request_samples), std::end(request_samples));
    std::vector<S> result {};
    result.reserve(request_samples.size());
    std::set_intersection(std::make_move_iterator(std::begin(request_samples)),
                          std::make_move_iterator(std::end(request_samples)),
                          std::cbegin(file_samples), std::cend(file_samples),
                          std::back_inserter(result));
    return result;
}

bool contains(const std::vector<HtslibSamFacade::SampleName>& samples, const HtslibSamFacade::SampleName& sample) noexcept
{
    return std::find(std::cbegin(samples), std::cend(samples), sample) != std::cend(samples);
}

} // namespace

// iterate

bool
HtslibSamFacade::iterate(const GenomicRegion& region,
                         AlignedReadReadVisitor visitor) const
{
    HtslibIterator itr {*this, region};
    if (samples_.size() == 1) {
        while (++itr) {
            if (!visitor(samples_.front(), *itr)) return false;
        }
    } else {
        while (++itr) {
            if (!visitor(sample_names_.at(itr.read_group()), *itr)) return false;
        }
    }
    return true;
}

bool
HtslibSamFacade::iterate(const SampleName& sample,
                         const GenomicRegion& region,
                         AlignedReadReadVisitor visitor) const
{
    if (contains(samples_, sample)) {
        if (samples_.size() == 1) {
            return iterate(region, visitor);
        } else {
            HtslibIterator itr {*this, region};
            while (++itr) {
                const auto& read_sample = sample_names_.at(itr.read_group());
                if (read_sample == sample) {
                    if (!visitor(sample, *itr)) return false;
                }
            }
        }
    }
    return true;
}

bool
HtslibSamFacade::iterate(const std::vector<SampleName>& samples,
                         const GenomicRegion& region,
                         AlignedReadReadVisitor visitor) const
{
    if (samples.empty()) return true;
    if (samples.size() == 1) {
        return iterate(samples.front(), region, visitor);
    } else if (set_equal(samples, samples_)) {
        return iterate(region, visitor);
    } else {
        const auto readable_samples = get_readable_samples(samples, samples_);
        HtslibIterator itr {*this, region};
        while (++itr) {
            const auto& sample = sample_names_.at(itr.read_group());
            if (std::binary_search(std::cbegin(readable_samples), std::cend(readable_samples), sample)) {
                if (!visitor(sample, *itr)) return false;
            }
        }
        return true;
    }
}

bool
HtslibSamFacade::iterate(const GenomicRegion& region,
                         ContigRegionVisitor visitor) const
{
    HtslibIterator itr {*this, region};
    if (samples_.size() == 1) {
        while (++itr) {
            if (!visitor(samples_.front(), itr.region())) return false;
        }
    } else {
        while (++itr) {
            if (!visitor(sample_names_.at(itr.read_group()), itr.region())) return false;
        }
    }
    return true;
}

bool
HtslibSamFacade::iterate(const SampleName& sample,
                         const GenomicRegion& region,
                         ContigRegionVisitor visitor) const
{
    if (contains(samples_, sample)) {
        if (samples_.size() == 1) {
            return iterate(region, visitor);
        } else {
            HtslibIterator itr {*this, region};
            while (++itr) {
                const auto& read_sample = sample_names_.at(itr.read_group());
                if (read_sample == sample) {
                    if (!visitor(sample, itr.region())) return false;
                }
            }
        }
    }
    return true;
}

bool
HtslibSamFacade::iterate(const std::vector<SampleName>& samples,
                         const GenomicRegion& region,
                         ContigRegionVisitor visitor) const
{
    if (samples.empty()) return true;
    if (samples.size() == 1) {
        return iterate(samples.front(), region, visitor);
    } else if (set_equal(samples, samples_)) {
        return iterate(region, visitor);
    } else {
        const auto readable_samples = get_readable_samples(samples, samples_);
        HtslibIterator itr {*this, region};
        while (++itr) {
            const auto& sample = sample_names_.at(itr.read_group());
            if (std::binary_search(std::cbegin(readable_samples), std::cend(readable_samples), sample)) {
                if (!visitor(sample, itr.region())) return false;
            }
        }
        return true;
    }
}

// has_reads

bool HtslibSamFacade::has_reads(const GenomicRegion& region) const
{
    HtslibIterator it {*this, region};
    return ++it;
}

bool HtslibSamFacade::has_reads(const SampleName& sample, const GenomicRegion& region) const
{
    if (samples_.size() == 1 && samples_.front() == sample) return has_reads(region);
    HtslibIterator it {*this, region};
    while (++it) if (sample_names_.at(it.read_group()) == sample) return true;
    return false;
}

bool HtslibSamFacade::has_reads(const std::vector<SampleName>& samples, const GenomicRegion& region) const
{
    if (samples.empty()) return false;
    if (samples.size() == 1) return has_reads(samples.front(), region);
    if (set_equal(samples, samples_)) return has_reads(region);
    const auto readable_samples = get_readable_samples(samples, samples_);
    HtslibIterator it {*this, region};
    while (++it) {
        if (std::binary_search(std::cbegin(readable_samples), std::cend(readable_samples),
                               sample_names_.at(it.read_group()))) {
            return true;
        }
    }
    return false;
}

// count_reads

std::size_t HtslibSamFacade::count_reads(const GenomicRegion& region) const
{
    HtslibIterator it {*this, region};
    std::size_t result {0};
    while (++it) ++result;
    return result;
}

std::size_t HtslibSamFacade::count_reads(const SampleName& sample, const GenomicRegion& region) const
{
    if (samples_.size() == 1 && samples_.front() == sample) return count_reads(region);
    HtslibIterator it {*this, region};
    std::size_t result {0};
    while (++it && sample_names_.at(it.read_group()) == sample) ++result;
    return result;
}

std::size_t HtslibSamFacade::count_reads(const std::vector<SampleName>& samples, const GenomicRegion& region) const
{
    if (samples.empty()) return 0;
    if (samples.size() == 1) return count_reads(samples.front(), region);
    if (set_equal(samples, samples_)) return count_reads(region);
    const auto readable_samples = get_readable_samples(samples, samples_);
    HtslibIterator it {*this, region};
    std::size_t result {0};
    while (++it && std::binary_search(std::cbegin(readable_samples), std::cend(readable_samples),
                                      sample_names_.at(it.read_group()))) {
        ++result;
    }
    return result;
}

// extract_read_positions

HtslibSamFacade::PositionList
HtslibSamFacade::extract_read_positions(const GenomicRegion& region, std::size_t max_coverage) const
{
    PositionList result {};
    result.reserve(max_coverage);
    HtslibIterator it {*this, region};
    while (max_coverage > 0 && ++it) {
        result.push_back(it.begin());
        --max_coverage;
    }
    return result;
}

HtslibSamFacade::PositionList
HtslibSamFacade::extract_read_positions(const SampleName& sample, const GenomicRegion& region,
                                        std::size_t max_coverage) const
{
    if (!contains(samples_, sample)) return {};
    if (samples_.size() == 1) return extract_read_positions(region, max_coverage);
    PositionList result {};
    result.reserve(max_coverage);
    HtslibIterator it {*this, region};
    while (max_coverage > 0 && ++it) {
        if (sample == sample_names_.at(it.read_group())) {
            result.push_back(it.begin());
            --max_coverage;
        }
    }
    return result;
}

HtslibSamFacade::PositionList
HtslibSamFacade::extract_read_positions(const std::vector<SampleName>& samples, const GenomicRegion& region,
                                        std::size_t max_coverage) const
{
    if (samples.empty()) return {};
    if (samples.size() == 1) return extract_read_positions(samples.front(), region, max_coverage);
    if (set_equal(samples, samples_)) return extract_read_positions(region, max_coverage);
    PositionList result {};
    result.reserve(max_coverage);
    HtslibIterator it {*this, region};
    while (max_coverage > 0 && ++it) {
        if (contains(samples, sample_names_.at(it.read_group()))) {
            result.push_back(it.begin());
            --max_coverage;
        }
    }
    return result;
}

// fetch_reads

namespace {

template <typename Container>
bool try_reserve(Container& c, const std::size_t max, const std::size_t min)
{
    assert(max >= min);
    if (max == 0) return true;
    for (auto curr = max; curr >= min; curr /= 2) {
        try {
            c.reserve(curr);
            return true;
        } catch (const std::bad_alloc& e) {}
    }
    return false;
}

} // namespace

HtslibSamFacade::SampleReadMap HtslibSamFacade::fetch_reads(const GenomicRegion& region) const
{
    SampleReadMap result {samples_.size()};
    if (samples_.size() == 1) {
        return {{samples_.front(), fetch_reads(samples_.front(), region)}};
    }
    HtslibIterator it {*this, region};
    for (const auto& sample : samples_) {
        auto p = result.emplace(std::piecewise_construct, std::forward_as_tuple(sample), std::forward_as_tuple());
        try_reserve(p.first->second, defaultReserve_, defaultReserve_ / 10);
    }
    while (++it) {
        try {
            result.at(sample_names_.at(it.read_group())).emplace_back(*it);
        } catch (InvalidBamRecord& e) {
            // TODO: Just ignore? Could log or something.
            //std::clog << "Warning: " << e.what() << std::endl;
        } catch (...) {
            throw;
        }
    }
    return result;
}

HtslibSamFacade::ReadContainer HtslibSamFacade::fetch_reads(const SampleName& sample, const GenomicRegion& region) const
{
    if (!contains(samples_, sample)) return {};
    if (samples_.size() == 1) return fetch_all_reads(region);
    HtslibIterator it {*this, region};
    ReadContainer result {};
    try_reserve(result, defaultReserve_, defaultReserve_ / 10);
    while (++it) {
        if (sample_names_.at(it.read_group()) == sample) {
            try {
                result.emplace_back(*it);
            } catch (InvalidBamRecord& e) {
                // TODO
            } catch (...) {
                throw;
            }
        }
    }
    return result;
}

HtslibSamFacade::SampleReadMap HtslibSamFacade::fetch_reads(const std::vector<SampleName>& samples,
                                                            const GenomicRegion& region) const
{
    if (samples.size() == 1) {
        return {{samples.front(), fetch_reads(samples.front(), region)}};
    }
    if (set_equal(samples_, samples)) return fetch_reads(region);
    HtslibIterator it {*this, region};
    SampleReadMap result {samples.size()};
    for (const auto& sample : samples) {
        if (contains(samples_, sample)) {
            auto p = result.emplace(std::piecewise_construct,
                                    std::forward_as_tuple(sample),
                                    std::forward_as_tuple());
            try_reserve(p.first->second, defaultReserve_, defaultReserve_ / 10);
        }
    }
    if (result.empty()) return result; // no matching samples
    while (++it) {
        const auto& sample = sample_names_.at(it.read_group());
        if (result.count(sample) == 1) {
            try {
                result.at(sample_names_.at(it.read_group())).emplace_back(*it);
            } catch (InvalidBamRecord& e) {
                // TODO
            } catch (...) {
                throw;
            }
        }
    }
    return result;
}

auto to_hts_regions(const std::vector<GenomicRegion>& regions)
{
    std::vector<hts_pair_pos_t> result(regions.size());
    std::transform(std::cbegin(regions), std::cend(regions), std::begin(result),
                   [] (const auto& region) -> hts_pair_pos_t { return {region.begin(), region.end()}; });
    return result;
}

HtslibSamFacade::SampleReadMap HtslibSamFacade::fetch_reads(const std::vector<GenomicRegion>& regions) const
{
    if (regions.size() == 1) {
        return fetch_reads(regions.front());
    }
    SampleReadMap result {samples_.size()};
    if (samples_.size() == 1) {
        return {{samples_.front(), fetch_reads(samples_.front(), regions)}};
    }
    auto hts_regions = to_hts_regions(regions);
    std::vector<hts_reglist_t> hts_region_list {make_hts_region_list(regions.front().contig_name(), hts_regions)};
    HtslibIterator it {*this, hts_region_list};
    for (const auto& sample : samples_) {
        auto p = result.emplace(std::piecewise_construct, std::forward_as_tuple(sample), std::forward_as_tuple());
        try_reserve(p.first->second, defaultReserve_, defaultReserve_ / 10);
    }
    while (++it) {
        try {
            result.at(sample_names_.at(it.read_group())).emplace_back(*it);
        } catch (InvalidBamRecord& e) {
            // TODO: Just ignore? Could log or something.
            //std::clog << "Warning: " << e.what() << std::endl;
        } catch (...) {
            throw;
        }
    }
    return result;
}

HtslibSamFacade::ReadContainer HtslibSamFacade::fetch_reads(const SampleName& sample, const std::vector<GenomicRegion>& regions) const
{
    if (regions.size() == 1) {
        return fetch_reads(sample, regions.front());
    }
    if (!contains(samples_, sample)) return {};
    if (samples_.size() == 1) return fetch_all_reads(regions);
    auto hts_regions = to_hts_regions(regions);
    std::vector<hts_reglist_t> hts_region_list {make_hts_region_list(regions.front().contig_name(), hts_regions)};
    HtslibIterator it {*this, hts_region_list};
    ReadContainer result {};
    try_reserve(result, defaultReserve_, defaultReserve_ / 10);
    while (++it) {
        if (sample_names_.at(it.read_group()) == sample) {
            try {
                result.emplace_back(*it);
            } catch (InvalidBamRecord& e) {
                // TODO
            } catch (...) {
                throw;
            }
        }
    }
    return result;
}

HtslibSamFacade::SampleReadMap HtslibSamFacade::fetch_reads(const std::vector<SampleName>& samples,
                                                            const std::vector<GenomicRegion>& regions) const
{
    if (regions.size() == 1) {
        return fetch_reads(samples, regions.front());
    }
    if (samples.size() == 1) {
        return {{samples.front(), fetch_reads(samples.front(), regions)}};
    }
    if (set_equal(samples_, samples)) return fetch_reads(regions);
    auto hts_regions = to_hts_regions(regions);
    std::vector<hts_reglist_t> hts_region_list {make_hts_region_list(regions.front().contig_name(), hts_regions)};
    HtslibIterator it {*this, hts_region_list};
    SampleReadMap result {samples.size()};
    for (const auto& sample : samples) {
        if (contains(samples_, sample)) {
            auto p = result.emplace(std::piecewise_construct,
                                    std::forward_as_tuple(sample),
                                    std::forward_as_tuple());
            try_reserve(p.first->second, defaultReserve_, defaultReserve_ / 10);
        }
    }
    if (result.empty()) return result; // no matching samples
    while (++it) {
        const auto& sample = sample_names_.at(it.read_group());
        if (result.count(sample) == 1) {
            try {
                result.at(sample_names_.at(it.read_group())).emplace_back(*it);
            } catch (InvalidBamRecord& e) {
                // TODO
            } catch (...) {
                throw;
            }
        }
    }
    return result;
}

std::vector<GenomicRegion::ContigName> HtslibSamFacade::reference_contigs() const
{
    std::vector<GenomicRegion::ContigName> result {};
    result.reserve(hts_header_->n_targets);
    for (HtsTid target {0}; target < hts_header_->n_targets; ++target) {
        result.emplace_back(get_contig_name(target));
    }
    return result;
}

boost::optional<std::vector<GenomicRegion::ContigName>> HtslibSamFacade::mapped_contigs() const
{
    std::vector<GenomicRegion::ContigName > result {};
    result.reserve(hts_header_->n_targets);
    for (HtsTid target {0}; target < hts_header_->n_targets; ++target) {
        const auto& target_contig_name = get_contig_name(target);
        // CRAM files don't seem to have the same index stats as BAM files so
        // we don't know which contigs have been mapped to
        if (hts_file_->is_cram || get_num_mapped_reads(target_contig_name) > 0) {
            result.push_back(target_contig_name);
        }
    }
    return result;
}

void HtslibSamFacade::write(const AlignedRead& read)
{
    if (!hts_file_ || !hts_header_) {
        throw UnwritableBAM {file_path_};
    }
    std::unique_ptr<bam1_t, HtsBam1Deleter> record {bam_init1(), HtsBam1Deleter {}};
    if (!record) {
        throw UnwritableBAM {file_path_};
    }
    write(read, record.get());
    if (sam_write1(hts_file_.get(), hts_header_.get(), record.get()) < 0) {
        throw UnwritableBAM {file_path_};
    }
}

// private methods

HtslibSamFacade::ReadContainer HtslibSamFacade::fetch_all_reads(const GenomicRegion& region) const
{
    HtslibIterator it {*this, region};
    ReadContainer result {};
    try_reserve(result, defaultReserve_, defaultReserve_ / 10);
    while (++it) {
        try {
            result.emplace_back(*it);
        } catch (InvalidBamRecord& e) {
            // TODO
        } catch (...) {
            throw;
        }
    }
    return result;
}

HtslibSamFacade::ReadContainer HtslibSamFacade::fetch_all_reads(const std::vector<GenomicRegion>& regions) const
{
    auto hts_regions = to_hts_regions(regions);
    std::vector<hts_reglist_t> hts_region_list {make_hts_region_list(regions.front().contig_name(), hts_regions)};
    HtslibIterator it {*this, hts_region_list};
    ReadContainer result {};
    try_reserve(result, defaultReserve_, defaultReserve_ / 10);
    while (++it) {
        try {
            result.emplace_back(*it);
        } catch (InvalidBamRecord& e) {
            // TODO
        } catch (...) {
            throw;
        }
    }
    return result;
}

hts_reglist_t HtslibSamFacade::make_hts_region_list(const GenomicRegion::ContigName& contig, std::vector<hts_pair_pos_t>& regions) const
{
    hts_reglist_t result {};
    result.reg = contig.c_str();
    result.tid = hts_targets_.at(contig);
    result.count = regions.size();
    result.intervals = regions.data();
    result.min_beg = regions.front().beg;
    result.max_end = regions.back().end;
    return result;
}

bool is_tag_type(const std::string& header_line, const std::string& tag)
{
    return header_line.compare(1, 2, tag) == 0;
}

bool has_tag(const std::string& header_line, const std::string& tag)
{
    return header_line.find(tag) != std::string::npos;
}

std::string extract_tag_value(const std::string& line, const std::string& tag)
{
    const auto tag_position = line.find(tag); // format is TAG:VALUE\t
    assert(tag_position != std::string::npos);
    const auto value_position = line.find(':', tag_position) + 1;
    const auto tag_value_size = line.find('\t', value_position) - value_position;
    return line.substr(value_position, tag_value_size);
}

void HtslibSamFacade::init_maps()
{
    hts_targets_.reserve(hts_header_->n_targets);
    contig_names_.reserve(hts_header_->n_targets);
    
    for (HtsTid target {0}; target < hts_header_->n_targets; ++target) {
        hts_targets_.emplace(hts_header_->target_name[target], target);
        contig_names_.emplace(target, hts_header_->target_name[target]);
    }
    
    const std::string header_text(hts_header_->text, hts_header_->l_text);
    std::istringstream header_ss {header_text};
    std::string line;
    unsigned num_read_groups {0};
    
    while (std::getline(header_ss, line, '\n')) {
        if (is_tag_type(line, readGroupTag)) {
            if (!has_tag(line, readGroupIdTag)) {
                MalformedBAMHeader e {file_path_};
                e.set_reason("a read group identifier tag (ID) in @RG lines is required but was not found");
                throw e;
            }
            if (!has_tag(line, sampleIdTag)) {
                // The SAM specification does not specify the sample tag 'SM' as a required,
                // however we can't do much without it.
                MalformedBAMHeader e {file_path_};
                e.set_reason("a sample tag (SM) in @RG lines is required but was not found");
                throw e;
            }
            sample_names_.emplace(extract_tag_value(line, readGroupIdTag),
                                  extract_tag_value(line, sampleIdTag));
            ++num_read_groups;
        }
    }
    
    if (num_read_groups == 0) {
        MalformedBAMHeader e {file_path_};
        e.set_reason("at least one read group (@RG) line is required but none were found");
        throw e;
    }
}

HtslibSamFacade::HtsTid HtslibSamFacade::get_htslib_target(const GenomicRegion::ContigName& contig) const
{
    return hts_targets_.at(contig);
}

const std::string& HtslibSamFacade::get_contig_name(HtsTid target) const
{
    return contig_names_.at(target);
}

// HtslibIterator

auto make_hts_iterator(const hts_idx_t* idx, bam_hdr_t* hdr, const GenomicRegion& region)
{
    const auto region_str = to_string(region);
    return sam_itr_querys(idx, hdr, region_str.c_str());
}

HtslibSamFacade::HtslibIterator::HtslibIterator(const HtslibSamFacade& hts_facade, const GenomicRegion& region)
: hts_facade_ {hts_facade}
, hts_iterator_ {hts_facade.is_open()
        ? make_hts_iterator(hts_facade_.hts_index_.get(), hts_facade_.hts_header_.get(), region)
    : nullptr, HtsIteratorDeleter {}}
, hts_bam1_ {bam_init1(), HtsBam1Deleter {}}
{
    if (hts_iterator_ == nullptr) {
        throw std::runtime_error {"HtslibIterator: could not load iterator for " + hts_facade.file_path_.string()};
    }
    if (hts_bam1_ == nullptr) {
        throw std::runtime_error {"HtslibIterator: error creating bam1 for " + hts_facade.file_path_.string()};
    }
}

HtslibSamFacade::HtslibIterator::HtslibIterator(const HtslibSamFacade& hts_facade, const GenomicRegion::ContigName& contig)
: hts_facade_ {hts_facade}
, hts_iterator_ {hts_facade.is_open() ? sam_itr_querys(hts_facade_.hts_index_.get(), hts_facade_.hts_header_.get(),
                                                     contig.c_str()) : nullptr, HtsIteratorDeleter {}}
, hts_bam1_ {bam_init1(), HtsBam1Deleter {}}
{
    if (hts_iterator_ == nullptr) {
        throw std::runtime_error {"HtslibIterator: could not load iterator for " + hts_facade.file_path_.string()};
    }
    if (hts_bam1_ == nullptr) {
        throw std::runtime_error {"HtslibIterator: error creating bam1 for " + hts_facade.file_path_.string()};
    }
}

HtslibSamFacade::HtslibIterator::HtslibIterator(const HtslibSamFacade& hts_facade, std::vector<hts_reglist_t>& regions)
: hts_facade_ {hts_facade}
, hts_iterator_ {hts_facade.is_open() ? sam_itr_regions(hts_facade_.hts_index_.get(), hts_facade_.hts_header_.get(),
                                                        regions.data(), regions.size()) : nullptr, HtsIteratorDeleter {}}
, hts_bam1_ {bam_init1(), HtsBam1Deleter {}}
{
    if (hts_iterator_ == nullptr) {
        throw std::runtime_error {"HtslibIterator: could not load iterator for " + hts_facade.file_path_.string()};
    }
    if (hts_bam1_ == nullptr) {
        throw std::runtime_error {"HtslibIterator: error creating bam1 for " + hts_facade.file_path_.string()};
    }
}

std::string extract_read_name(const bam1_t* b)
{
    return std::string {bam_get_qname(b)};
}

bool HtslibSamFacade::HtslibIterator::operator++()
{
    return sam_itr_next(hts_facade_.hts_file_.get(), hts_iterator_.get(), hts_bam1_.get()) >= 0;
}

auto extract_read_pos(const bam1_t* b) noexcept
{
    return b->core.pos;
}

auto extract_sequence_length(const bam1_t* b) noexcept
{
    return b->core.l_qseq;
}

char extract_base(const std::uint8_t* hts_sequence, const std::uint32_t index) noexcept
{
    static constexpr const char* symbolTable {"=ACMGRSVTWYHKDBN"};
    return symbolTable[bam_seqi(hts_sequence, index)];
}

AlignedRead::NucleotideSequence extract_sequence(const bam1_t* b)
{
    using NucleotideSequence = AlignedRead::NucleotideSequence;
    const auto sequence_length  = static_cast<NucleotideSequence::size_type>(extract_sequence_length(b));
    const auto hts_sequence     = bam_get_seq(b);
    NucleotideSequence result(sequence_length, 'N');
    std::uint32_t i {0};
    std::generate_n(std::begin(result), sequence_length, [&] () { return extract_base(hts_sequence, i++); });
    return result;
}

AlignedRead::BaseQualityVector extract_qualities(const bam1_t* b)
{
    const auto qualities = bam_get_qual(b);
    const auto length    = extract_sequence_length(b);
    return AlignedRead::BaseQualityVector(qualities, qualities + length);
}

auto get_cigar_length(const bam1_t* b) noexcept
{
    return b->core.n_cigar;
}

CigarString extract_cigar_string(const bam1_t* b)
{
    const auto cigar_operations = bam_get_cigar(b);
    const auto cigar_length     = get_cigar_length(b);
    CigarString result(cigar_length);
    std::transform(cigar_operations, cigar_operations + cigar_length, std::begin(result),
                   [] (const auto op) noexcept {
                       return CigarOperation {static_cast<CigarOperation::Size>(bam_cigar_oplen(op)),
                                              static_cast<CigarOperation::Flag>(bam_cigar_opchr(op))};
                   });
    return result;
}

// Some of these flags will need to be changes when htslib catches up to the new SAM spec
AlignedRead::Flags extract_flags(const bam1_core_t& c) noexcept
{
    AlignedRead::Flags result {};
    result.multiple_segment_template    = (c.flag & BAM_FPAIRED)        != 0;
    result.all_segments_in_read_aligned = (c.flag & BAM_FPROPER_PAIR)   != 0;
    result.unmapped                     = (c.flag & BAM_FUNMAP)         != 0;
    result.reverse_mapped               = (c.flag & BAM_FREVERSE)       != 0;
    result.secondary_alignment          = (c.flag & BAM_FSECONDARY)     != 0;
    result.qc_fail                      = (c.flag & BAM_FQCFAIL)        != 0;
    result.duplicate                    = (c.flag & BAM_FDUP)           != 0;
    result.supplementary_alignment      = (c.flag & BAM_FSUPPLEMENTARY) != 0;
    result.first_template_segment       = (c.flag & BAM_FREAD1)         != 0;
    result.last_template_segment        = (c.flag & BAM_FREAD2)         != 0;
    return result;
}

auto mapping_quality(const bam1_core_t& c) noexcept
{
    return static_cast<AlignedRead::MappingQuality>(c.qual);
}

bool has_multiple_segments(const bam1_core_t& c) noexcept
{
    return c.mtid != -1;
}

auto next_segment_position(const bam1_core_t& c) noexcept
{
    return static_cast<AlignedRead::MappingDomain::Position>(c.mpos);
}

auto template_length(const bam1_core_t& c) noexcept
{
    return static_cast<AlignedRead::MappingDomain::Size>(std::abs(c.isize));
}

auto extract_next_segment_flags(const bam1_core_t& c) noexcept
{
    AlignedRead::Segment::Flags result {};
    result.unmapped       = (c.flag & BAM_FMUNMAP)   != 0;
    result.reverse_mapped = (c.flag & BAM_FMREVERSE) != 0;
    return result;
}

auto get_annotations(const bam1_t* record, const std::vector<AlignedRead::Tag>& tags)
{
    std::vector<std::pair<AlignedRead::Tag, AlignedRead::Annotation>> result {};
    result.reserve(tags.size());
    for (const auto& tag : tags) {
        const auto ptr = bam_aux_get(record, tag.data());
        if (ptr) result.emplace_back(tag, bam_aux2Z(ptr));
    }
    return result;
}

static inline int aux_type2size(const std::uint8_t type)
{
    switch (type) {
    case 'A': case 'c': case 'C':
        return 1;
    case 's': case 'S':
        return 2;
    case 'i': case 'I': case 'f':
        return 4;
    case 'd':
        return 8;
    case 'Z': case 'H': case 'B':
        return type;
    default:
        return 0;
    }
}

auto get_annotations(const bam1_t* record)
{
    auto s = bam_get_aux(record);
    const auto end = record->data + record->l_data;
    std::vector<std::pair<AlignedRead::Tag, AlignedRead::Annotation>> result {};
    AlignedRead::Tag tag {};
    AlignedRead::Annotation annotation {};
    while (s != nullptr && end - s >= 3) {
        std::copy_n(s, tag.size(), std::begin(tag));
        s += tag.size();
        const auto tagtype = static_cast<char>(*s);
        switch (tagtype) {
            case 'c':
                // fall through
            case 'C':
                annotation = std::to_string(static_cast<int>(bam_aux2i(s)));
                s += 1;
                break;
            case 's':
                // fall through
            case 'S':
                annotation = std::to_string(static_cast<int>(bam_aux2i(s)));
                s += 2;
                break;
            case 'i':
                // fall through
            case 'I':
                annotation = std::to_string(static_cast<int>(bam_aux2i(s)));
                s += 4;
                break;
            case 'f':
                annotation = std::to_string(static_cast<float>(bam_aux2f(s)));
                s += 4;
                break;
            case 'd':
                annotation = std::to_string(static_cast<double>(bam_aux2f(s)));
                s += 8;
                break;
            case 'a':
                // fall through
            case 'A':
                annotation = std::to_string(static_cast<char>(bam_aux2A(s)));
                s += 1;
                break;
            case 'Z':
                // fall through
            case 'H':
                annotation = std::string {bam_aux2Z(s)};
                s += (annotation.length() + 1);
                break;
            case 'B': {
                const auto size = aux_type2size(*s);
                ++s;
                const auto n = le_to_u32(s);
                s += 4;
                s += (n * size);
                break;
            }
            default:
                throw std::runtime_error {"Unknown BAM tag type"};
        }
        result.emplace_back(tag, annotation);
        s += 1;
    }
    return result;
}

AlignedRead HtslibSamFacade::HtslibIterator::operator*() const
{
    using std::begin; using std::end; using std::next; using std::move;
    auto qualities = extract_qualities(hts_bam1_.get());
    auto cigar = extract_cigar_string(hts_bam1_.get());
    const auto& info = hts_bam1_->core;
    auto read_begin_tmp = clipped_begin(cigar, info.pos);
    auto sequence = extract_sequence(hts_bam1_.get());
    if (sequence.size() != qualities.size()) {
        throw InvalidBamRecord {hts_facade_.file_path_, extract_read_name(hts_bam1_.get()), "corrupt sequence data"};
    }
    if (read_begin_tmp < 0) {
        // Then the read hangs off the left of the contig, and we must remove bases, base_qualities, and
        // adjust the cigar string as we cannot have a negative begin position
        const auto overhang_size = static_cast<unsigned>(std::abs(read_begin_tmp));
        if (!qualities.empty() && qualities.size() == sequence.size()) {
            sequence.erase(begin(sequence), next(begin(sequence), overhang_size));
            qualities.erase(begin(qualities), next(begin(qualities), overhang_size));
        }
        auto soft_clip_size = cigar.front().size();
        if (overhang_size == soft_clip_size) {
            cigar.erase(begin(cigar));
        } else { // then soft_clip_size > overhang_size
            cigar.front() = CigarOperation {soft_clip_size - overhang_size, CigarOperation::Flag::softClipped};
        }
        read_begin_tmp = 0;
    }
    const auto read_begin = static_cast<AlignedRead::MappingDomain::Position>(read_begin_tmp);
    const auto& contig_name = hts_facade_.get_contig_name(info.tid);
    AlignedRead result;
    if (has_multiple_segments(info)) {
        result = {
            extract_read_name(hts_bam1_.get()),
            GenomicRegion {contig_name, read_begin, read_begin + octopus::reference_size<AlignedRead::MappingDomain::Position>(cigar)},
            move(sequence),
            move(qualities),
            move(cigar),
            mapping_quality(info),
            extract_flags(info),
            read_group(),
            hts_facade_.get_contig_name(info.mtid),
            next_segment_position(info),
            template_length(info),
            extract_next_segment_flags(info),
            get_annotations(hts_bam1_.get())
        };
    } else {
        result = {
            extract_read_name(hts_bam1_.get()),
            GenomicRegion {contig_name, read_begin, read_begin + octopus::reference_size<AlignedRead::MappingDomain::Size>(cigar)},
            move(sequence),
            move(qualities),
            move(cigar),
            mapping_quality(info),
            extract_flags(info),
            read_group(),
            get_annotations(hts_bam1_.get())
        };
    }
    return result;
}

HtslibSamFacade::ReadGroupIdType HtslibSamFacade::HtslibIterator::read_group() const
{
    const auto ptr = bam_aux_get(hts_bam1_.get(), readGroupTag.c_str());
    if (ptr == nullptr) {
        throw InvalidBamRecord {hts_facade_.file_path_, extract_read_name(hts_bam1_.get()), "no read group"};
    }
    return HtslibSamFacade::ReadGroupIdType {bam_aux2Z(ptr)};
}

bool HtslibSamFacade::HtslibIterator::is_good() const noexcept
{
    if (extract_sequence_length(hts_bam1_.get()) == 0) {
        return false;
    }
    if (bam_get_qual(hts_bam1_.get())[0] == 0xff) {
        return false;
    }
    const auto cigar_length = get_cigar_length(hts_bam1_.get());
    if (cigar_length == 0) return false;
    const auto cigar_operations = bam_get_cigar(hts_bam1_.get());
    return std::all_of(cigar_operations, cigar_operations + cigar_length,
                       [] (const auto op) { return bam_cigar_oplen(op) > 0; });
}

ContigRegion HtslibSamFacade::HtslibIterator::region() const
{
    const auto cigar = extract_cigar_string(hts_bam1_.get());
    const auto begin = clipped_begin(cigar, static_cast<ContigRegion::Position>(hts_bam1_->core.pos));
    return ContigRegion {begin, begin + octopus::reference_size(cigar)};
}

std::size_t HtslibSamFacade::HtslibIterator::begin() const noexcept
{
    return hts_bam1_->core.pos;
}

namespace {

void set_contig(const std::int32_t tid, bam1_t* result) noexcept
{
    result->core.tid = tid;
}

void set_pos(const AlignedRead& read, bam1_t* result) noexcept
{
    result->core.pos = mapped_begin(read);
    if (is_front_soft_clipped(read)) {
        result->core.pos += get_soft_clipped_sizes(read).first;
    }
}

void set_mapping_quality(const AlignedRead& read, bam1_t* result) noexcept
{
    result->core.qual = read.mapping_quality();
}

void set_flag(bool set, std::uint16_t mask, std::uint16_t& result) noexcept
{
    constexpr std::uint16_t zeros {0}, ones {std::numeric_limits<std::uint16_t>::max()};
    result |= (set ? ones : zeros) & mask;
}

void set_flags(const AlignedRead& read, bam1_t* result) noexcept
{
    auto& bitset = result->core.flag;
    set_flag(read.is_marked_multiple_segment_template(),    BAM_FPAIRED, bitset);
    set_flag(read.is_marked_all_segments_in_read_aligned(), BAM_FPROPER_PAIR, bitset);
    set_flag(read.is_marked_unmapped(),                     BAM_FUNMAP, bitset);
    set_flag(read.is_marked_next_segment_unmapped(),        BAM_FMUNMAP, bitset);
    set_flag(read.is_marked_reverse_mapped(),               BAM_FREVERSE, bitset);
    set_flag(read.is_marked_next_segment_reverse_mapped(),  BAM_FMREVERSE, bitset);
    set_flag(read.is_marked_secondary_alignment(),          BAM_FSECONDARY, bitset);
    set_flag(read.is_marked_qc_fail(),                      BAM_FQCFAIL, bitset);
    set_flag(read.is_marked_duplicate(),                    BAM_FDUP, bitset);
    set_flag(read.is_marked_supplementary_alignment(),      BAM_FSUPPLEMENTARY, bitset);
    set_flag(read.is_marked_first_template_segment(),       BAM_FREAD1, bitset);
    set_flag(read.is_marked_last_template_segment(),        BAM_FREAD2, bitset);
}

void set_segment(const AlignedRead& read, const std::int32_t tid, bam1_t* result)
{
    const auto& segment = read.next_segment();
    result->core.mtid = tid;
    result->core.mpos = segment.begin();
    result->core.isize = static_cast<std::int32_t>(segment.inferred_template_length());
    if (segment.begin() < mapped_begin(read)) {
        result->core.isize *= -1;
    }
    set_flag(segment.is_marked_unmapped(),       BAM_FMUNMAP, result->core.flag);
    set_flag(segment.is_marked_reverse_mapped(), BAM_FMREVERSE, result->core.flag);
}

auto name_bytes(const AlignedRead& read) noexcept
{
    auto result = read.name().size() + 1;
    if (result % 4 > 0) result +=  4 - (result % 4); // next address must be 4-byte aligned
    return result;
}

auto sequence_bytes(const AlignedRead& read) noexcept
{
    return (read.sequence().size() + 1) / 2; // 4 bits per base
}

auto base_quality_bytes(const AlignedRead& read) noexcept
{
    return read.base_qualities().size();
}

auto cigar_bytes(const AlignedRead& read) noexcept
{
    return 4 * read.cigar().size();
}

auto calculate_required_field_bytes(const AlignedRead& read) noexcept
{
    std::size_t result {0};
    result += name_bytes(read);
    result += sequence_bytes(read);
    result += base_quality_bytes(read);
    result += cigar_bytes(read);
    return result;
}


auto calculate_annotation_bytes(const AlignedRead::Tag& tag, const AlignedRead::Annotation& annotation)
{
    return tag.size() + annotation.size() + 2;// 1 for tag, 1 for '\0'
}

auto calculate_aux_bytes(const AlignedRead& read)
{
    std::size_t result {0};
    for (auto tag : read.tags()) {
        result += calculate_annotation_bytes(tag, *read.annotation(tag));
    }
    return result;
}

template <typename Read>
void allocate_variable_length_data(const Read& read, bam1_t* result)
{
    result->m_data = calculate_required_field_bytes(read) + calculate_aux_bytes(read);
    result->data = (std::uint8_t*) std::realloc(result->data, result->m_data);
    std::fill_n(result->data, result->m_data, 0);
}

void set_name(const AlignedRead& read, bam1_t* result)
{
    const auto& name = read.name();
    std::copy(std::cbegin(name), std::cend(name), result->data);
    result->core.l_qname = name.size() + 1;
    result->core.l_extranul = (result->core.l_qname % 4 > 0) ? 4 - (result->core.l_qname % 4) : 0; // next address must be 4 byte aligned
    result->core.l_qname += result->core.l_extranul;
    result->l_data += result->core.l_qname;
}

void set_cigar(const AlignedRead& read, bam1_t* result) noexcept
{
    const auto& cigar = read.cigar();
    result->core.n_cigar = cigar.size();
    std::transform(std::cbegin(cigar), std::cend(cigar), bam_get_cigar(result),
                   [] (const CigarOperation& op) noexcept -> std::uint32_t {
                       // Lower 4 bits for CIGAR operation and the higher 28 bits for size
                       std::uint32_t result = op.size();
                       result <<= BAM_CIGAR_SHIFT;
                       using Flag = CigarOperation::Flag;
                       switch (op.flag()) {
                           case Flag::alignmentMatch: result |= BAM_CMATCH; break;
                           case Flag::insertion:      result |= BAM_CINS; break;
                           case Flag::deletion:       result |= BAM_CDEL; break;
                           case Flag::skipped:        result |= BAM_CREF_SKIP; break;
                           case Flag::softClipped:    result |= BAM_CSOFT_CLIP; break;
                           case Flag::hardClipped:    result |= BAM_CHARD_CLIP; break;
                           case Flag::padding:        result |= BAM_CPAD; break;
                           case Flag::sequenceMatch:  result |= BAM_CEQUAL; break;
                           case Flag::substitution:   result |= BAM_CDIFF; break;
                       }
                       return result;
                   });
    result->l_data += cigar_bytes(read);
}

static constexpr std::array<std::uint8_t, 128> sam_bases
{
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 1, 0, 2, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 15,0,
0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

void set_read_sequence(const AlignedRead& read, bam1_t* result) noexcept
{
    const auto& sequence = read.sequence();
    result->core.l_qseq = sequence.size();
    // Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
    // 8 for T and 15 for N. Two bases are packed in one byte with the base
    // at the higher 4 bits having smaller coordinate on the read.
    auto bam_seq_itr = bam_get_seq(result);
    for (std::size_t i {1}; i < sequence.size(); i += 2, ++bam_seq_itr) {
        *bam_seq_itr = sam_bases[sequence[i - 1]] << 4 | sam_bases[sequence[i]];
    }
    if (sequence.size() % 2 == 1) {
        *bam_seq_itr = sam_bases[sequence.back()] << 4;
    }
    result->l_data += sequence_bytes(read);
}

void set_base_qualities(const AlignedRead& read, bam1_t* result) noexcept
{
    std::copy(std::cbegin(read.base_qualities()), std::cend(read.base_qualities()), bam_get_qual(result));
    result->l_data += base_quality_bytes(read);
}

void set_required_data(const AlignedRead& read, bam1_t* result) noexcept
{
    set_name(read, result);
    set_cigar(read, result);
    set_read_sequence(read, result);
    set_base_qualities(read, result);
}

void set_annotation(const AlignedRead::Tag& tag, const AlignedRead::Annotation& annotation, bam1_t* result) noexcept
{
    // TODO - detect other data types from annotation string
    bam_aux_update_str(result, tag.data(), annotation.size() + 1, annotation.c_str());
}

void set_aux_data(const AlignedRead& read, bam1_t* result) noexcept
{
    for (auto tag : read.tags()) {
        set_annotation(tag, *read.annotation(tag), result);
    }
}

template <typename Read>
void set_variable_length_data(const Read& read, bam1_t* result)
{
    allocate_variable_length_data(read, result);
    set_required_data(read, result);
    set_aux_data(read, result);
}

} // namespace

void HtslibSamFacade::set_fixed_length_data(const AlignedRead& read, bam1_t* result) const
{
    set_contig(hts_targets_.at(contig_name(read)), result);
    set_pos(read, result);
    set_mapping_quality(read, result);
    set_flags(read, result);
    if (read.has_other_segment()) {
        set_segment(read, hts_targets_.at(read.next_segment().contig_name()), result);
    } else {
        result->core.mtid = -1;
        result->core.mpos = 0;
        result->core.isize = 0;
    }
}

void HtslibSamFacade::write(const AlignedRead& read, bam1_t* result) const
{
    set_fixed_length_data(read, result);
    set_variable_length_data(read, result);
}

} // namespace io
} // namespace octopus
