// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "htslib_sam_facade.hpp"

#include <sstream>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <sstream>

#include <boost/filesystem/operations.hpp>

#include <basics/cigar_string.hpp>
#include <basics/genomic_region.hpp>
#include <basics/contig_region.hpp>
#include <exceptions/missing_file_error.hpp>
#include <exceptions/missing_index_error.hpp>
#include <exceptions/malformed_file_error.hpp>

namespace octopus { namespace io {

static const std::string readGroupTag   {"RG"};
static const std::string readGroupIdTag {"ID"};
static const std::string sampleIdTag    {"SM"};

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

class InvalidBamRecord : public std::runtime_error
{
public:
    InvalidBamRecord(boost::filesystem::path file_path, std::string read_name, std::string message)
    : runtime_error {"Invalid BAM record"}
    , file_path_ {std::move(file_path)}
    , read_name_ {std::move(read_name)}
    , message_ {std::move(message)}
    , msg_ {}
    {}
    
    virtual ~InvalidBamRecord() = default;
    
    virtual const char* what() const noexcept override
    {
        std::ostringstream ss {};
        ss << runtime_error::what() << ": in file " << file_path_ << ", in read " << read_name_
                << ": " << message_;
        msg_ = ss.str();
        return msg_.c_str();
    }
private:
    boost::filesystem::path file_path_;
    std::string message_, read_name_;
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
}

HtslibSamFacade::HtslibSamFacade(Path file_path)
: file_path_ {std::move(file_path)}
, hts_file_ {open_hts_file(file_path_), HtsFileDeleter {}}
, hts_header_ {(hts_file_) ? sam_hdr_read(hts_file_.get()) : nullptr, HtsHeaderDeleter {}}
, hts_index_ {(hts_file_) ? sam_index_load(hts_file_.get(), file_path_.c_str()) : nullptr, HtsIndexDeleter {}}
, hts_tids_ {}
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

unsigned HtslibSamFacade::count_reference_contigs() const
{
    if (!is_open()) {
        throw std::runtime_error {"HtslibSamFacade: attempting operation on closed file"};
    }
    return hts_header_->n_targets;
}

ContigRegion::Size HtslibSamFacade::get_reference_contig_size(const std::string& contig_name) const
{
    return hts_header_->target_len[get_htslib_tid(contig_name)];
}

std::uint64_t HtslibSamFacade::get_num_mapped_reads(const std::string& contig_name) const
{
    uint64_t num_mapped {}, num_unmapped {};
    hts_idx_get_stat(hts_index_.get(), get_htslib_tid(contig_name), &num_mapped, &num_unmapped);
    return num_mapped;
}

std::vector<HtslibSamFacade::SampleName> HtslibSamFacade::extract_samples() const
{
    return samples_;
}

std::vector<std::string> HtslibSamFacade::extract_read_groups_in_sample(const SampleName& sample) const
{
    std::vector<std::string> result {};
    
    for (const auto pair : sample_names_) {
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
bool is_subset_helper(const std::vector<T>& lhs, const std::vector<T>& rhs)
{
    return std::includes(cbegin(lhs), cend(lhs), cbegin(rhs), cend(rhs));
}

template <typename T>
bool is_subset_helper_sort_rhs(const std::vector<T>& lhs, std::vector<T> rhs)
{
    std::sort(std::begin(rhs), std::end(rhs));
    return is_subset_helper(lhs, rhs);
}

template <typename T>
bool is_subset_helper_sort_lhs(std::vector<T> lhs, const std::vector<T>& rhs)
{
    std::sort(std::begin(lhs), std::end(lhs));
    return is_subset_helper(lhs, rhs);
}

template <typename T>
bool is_subset_helper_sort_both(std::vector<T> lhs, std::vector<T> rhs)
{
    std::sort(std::begin(lhs), std::end(lhs));
    std::sort(std::begin(rhs), std::end(rhs));
    return is_subset_helper(lhs, rhs);
}

// Returns true iff rhs a subset of lhs.
template <typename T>
bool is_subset(const std::vector<T>& lhs, const std::vector<T>& rhs)
{
    if (is_sorted(lhs)) {
        return is_sorted(rhs) ? is_subset_helper(lhs, rhs) : is_subset_helper_sort_rhs(lhs, rhs);
    } else if (is_sorted(rhs)) {
        return is_subset_helper_sort_lhs(lhs, rhs);
    } else {
        return is_subset_helper_sort_both(lhs, rhs);
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

} // namespace

// has_reads

bool HtslibSamFacade::has_reads(const GenomicRegion& region) const
{
    HtslibIterator it {*this, region};
    return ++it;
}

bool HtslibSamFacade::has_reads(const SampleName& sample,
                                const GenomicRegion& region) const
{
    if (samples_.size() == 1 && samples_.front() == sample) return has_reads(region);
    HtslibIterator it {*this, region};
    while (++it) if (sample_names_.at(it.read_group()) == sample) return true;
    return false;
}

bool HtslibSamFacade::has_reads(const std::vector<SampleName>& samples,
                                const GenomicRegion& region) const
{
    if (samples.empty()) return false;
    if (samples.size() == 1) return has_reads(samples.front(), region);
    if (is_subset(samples, samples_)) return has_reads(region);
    
    const auto readable_samples = get_readable_samples(samples, samples_);
    HtslibIterator it {*this, region};
    
    while (++it) {
        if (std::binary_search(std::cbegin(readable_samples),
                               std::cend(readable_samples),
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

std::size_t HtslibSamFacade::count_reads(const SampleName& sample,
                                         const GenomicRegion& region) const
{
    if (samples_.size() == 1 && samples_.front() == sample) return count_reads(region);
    HtslibIterator it {*this, region};
    std::size_t result {0};
    while (++it && sample_names_.at(it.read_group()) == sample) ++result;
    return result;
}

std::size_t HtslibSamFacade::count_reads(const std::vector<SampleName>& samples,
                                         const GenomicRegion& region) const
{
    if (samples.empty()) return 0;
    if (samples.size() == 1) return count_reads(samples.front(), region);
    if (is_subset(samples, samples_)) return count_reads(region);
    
    const auto readable_samples = get_readable_samples(samples, samples_);
    HtslibIterator it {*this, region};
    std::size_t result {0};
    
    while (++it && std::binary_search(std::cbegin(readable_samples),
                                      std::cend(readable_samples),
                                      sample_names_.at(it.read_group()))) ++result;
    
    return result;
}
        
// find_covered_subregion

GenomicRegion expand_subregion(const GenomicRegion& region, const std::size_t remaining_coverage_quota,
                               const std::size_t first_read_begin, const std::size_t last_read_begin,
                               std::vector<unsigned>& position_counts)
{
    assert(first_read_begin <= last_read_begin);
    const auto result_begin = std::min(first_read_begin, static_cast<std::size_t>(region.begin()));
    const auto result_end   = (remaining_coverage_quota > 0) ? region.end() : last_read_begin;
    assert(result_end >= result_begin);
    position_counts.resize(result_end - result_begin, 0);
    
    if (result_begin < first_read_begin) {
        std::rotate(std::begin(position_counts),
                    std::next(std::begin(position_counts), last_read_begin - first_read_begin),
                    std::end(position_counts));
    }
    
    return GenomicRegion {
        region.contig_name(),
        static_cast<GenomicRegion::Size>(result_begin),
        static_cast<GenomicRegion::Size>(result_end)
    };
}

HtslibSamFacade::CoveragePair
HtslibSamFacade::find_covered_subregion(const GenomicRegion& region, std::size_t max_coverage) const
{
    HtslibIterator it {*this, region};
    
    if (max_coverage == 0 || !++it) {
        return std::make_pair(head_region(region), std::vector<unsigned> {});
    }
    
    --max_coverage;
    auto first_begin = it.begin();
    std::vector<unsigned> position_counts(max_coverage, 0);
    ++position_counts[0];
    auto read_begin = first_begin;
    
    while (++it && --max_coverage > 0) {
        read_begin = it.begin();
        const auto offset = read_begin - first_begin;
        if (offset >= position_counts.size()) {
            position_counts.resize(std::max(offset, position_counts.size() + max_coverage));
        } else {
            ++position_counts[offset];
        }
    }
    
    auto result_region = expand_subregion(region, max_coverage, first_begin, read_begin, position_counts);
    return std::make_pair(std::move(result_region), std::move(position_counts));
}

bool contains(const std::vector<HtslibSamFacade::SampleName>& samples,
              const HtslibSamFacade::SampleName& sample)
{
    return std::find(std::cbegin(samples), std::cend(samples), sample) != std::cend(samples);
}

HtslibSamFacade::CoveragePair
HtslibSamFacade::find_covered_subregion(const SampleName& sample,
                                        const GenomicRegion& region,
                                        std::size_t max_coverage) const
{
    if (!contains(samples_, sample)) {
        return std::make_pair(head_region(region), std::vector<unsigned> {});
    }
    
    HtslibIterator it {*this, region};
    
    if (max_coverage == 0 || !++it) {
        return std::make_pair(head_region(region), std::vector<unsigned> {});
    }
    if (samples_.size() > 1) {
        while (sample != sample_names_.at(it.read_group())) {
            if (!++it) return std::make_pair(head_region(region), std::vector<unsigned> {});
        };
    }
    
    --max_coverage;
    const auto first_begin = it.begin();
    std::vector<unsigned> position_counts(max_coverage, 0);
    ++position_counts[0];
    auto read_begin = first_begin;
    
    if (samples_.size() == 1) {
        while (++it && --max_coverage > 0) {
            read_begin = it.begin();
            const auto offset = read_begin - first_begin;
            if (offset >= position_counts.size()) {
                position_counts.resize(std::max(offset, position_counts.size() + max_coverage));
            } else {
                ++position_counts[offset];
            }
        }
    } else {
        while (++it && sample == sample_names_.at(it.read_group()) && --max_coverage > 0) {
            read_begin = it.begin();
            const auto offset = read_begin - first_begin;
            if (offset >= position_counts.size()) {
                position_counts.resize(std::max(offset, position_counts.size() + max_coverage));
            } else {
                ++position_counts[offset];
            }
        }
    }
    
    auto result_region = expand_subregion(region, max_coverage, first_begin, read_begin, position_counts);
    
    return std::make_pair(std::move(result_region), std::move(position_counts));
}

HtslibSamFacade::CoveragePair
HtslibSamFacade::find_covered_subregion(const std::vector<SampleName>& samples,
                                        const GenomicRegion& region,
                                        std::size_t max_coverage) const
{
    if (samples.empty()) return {head_region(region), {}};
    if (samples.size() == 1) {
        return find_covered_subregion(samples.front(), region, max_coverage);
    }
    if (is_subset(samples, samples_)) {
        return find_covered_subregion(region, max_coverage);
    }
    
    HtslibIterator it {*this, region};
    
    if (max_coverage == 0 || !++it) {
        return std::make_pair(head_region(region), std::vector<unsigned> {});
    }
    
    while (!contains(samples, sample_names_.at(it.read_group()))) {
        if (!++it) return std::make_pair(head_region(region), std::vector<unsigned> {});
    };
    
    --max_coverage;
    
    auto first_begin = it.begin();
    
    std::vector<unsigned> position_counts(max_coverage, 0);
    
    ++position_counts[0];
    
    auto read_begin = first_begin;
    
    while (++it && !contains(samples, sample_names_.at(it.read_group())) && --max_coverage > 0) {
        read_begin = it.begin();
        const auto offset = read_begin - first_begin;
        if (offset >= position_counts.size()) {
            position_counts.resize(std::max(offset, position_counts.size() + max_coverage));
        } else {
            ++position_counts[offset];
        }
    }
    
    auto result_region = expand_subregion(region, max_coverage, first_begin, read_begin, position_counts);
    return std::make_pair(std::move(result_region), std::move(position_counts));
}

// fetch_reads

HtslibSamFacade::SampleReadMap HtslibSamFacade::fetch_reads(const GenomicRegion& region) const
{
    SampleReadMap result {samples_.size()};
    
    if (samples_.size() == 1) {
        auto reads = fetch_reads(samples_.front(), region);
        result.emplace(samples_.front(), std::move(reads));
        return result;
    }
    
    HtslibIterator it {*this, region};
    
    for (const auto& sample : samples_) {
        auto p = result.emplace(std::piecewise_construct,
                                std::forward_as_tuple(sample),
                                std::forward_as_tuple());
        p.first->second.reserve(defaultReserve_);
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

HtslibSamFacade::ReadContainer HtslibSamFacade::fetch_reads(const SampleName& sample,
                                                            const GenomicRegion& region) const
{
    HtslibIterator it {*this, region};
    
    ReadContainer result {};
    
    if (std::find(std::cbegin(samples_), std::cend(samples_), sample) == std::cend(samples_)) {
        return result;
    }
    
    result.reserve(defaultReserve_);
    
    if (samples_.size() == 1) {
        while (++it) {
            try {
                result.emplace_back(*it);
            } catch (InvalidBamRecord& e) {
                // TODO
            } catch (...) {
                throw;
            }
        }
    } else {
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
    }
    
    return result;
}

HtslibSamFacade::SampleReadMap HtslibSamFacade::fetch_reads(const std::vector<SampleName>& samples,
                                                            const GenomicRegion& region) const
{
    if (is_subset(samples_, samples)) return fetch_reads(region);
    
    HtslibIterator it {*this, region};
    SampleReadMap result {samples.size()};
    
    for (const auto& sample : samples) {
        if (std::find(std::cbegin(samples_), std::cend(samples_), sample) != std::cend(samples_)) {
            auto p = result.emplace(std::piecewise_construct,
                                    std::forward_as_tuple(sample),
                                    std::forward_as_tuple());
            p.first->second.reserve(defaultReserve_);
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

std::vector<std::string> HtslibSamFacade::extract_reference_contig_names() const
{
    std::vector<std::string> result {};
    result.reserve(hts_header_->n_targets);
    
    for (HtsTidType hts_tid {}; hts_tid < hts_header_->n_targets; ++hts_tid) {
        result.emplace_back(contig_name(hts_tid));
    }
    
    return result;
}

std::vector<GenomicRegion> HtslibSamFacade::extract_possible_regions_in_file() const
{
    std::vector<GenomicRegion> result {};
    result.reserve(hts_header_->n_targets);
    
    for (HtsTidType hts_tid {0}; hts_tid < hts_header_->n_targets; ++hts_tid) {
        const auto& tid_contig_name = contig_name(hts_tid);
        // CRAM files don't seem to have the same index stats as BAM files so
        // we don't know which contigs have been mapped to
        if (hts_file_->is_cram || get_num_mapped_reads(tid_contig_name) > 0) {
            // without actually looking at the reads we don't know what coverage there is so
            // we just assume the entire contig is possible
            result.emplace_back(tid_contig_name, 0, get_reference_contig_size(tid_contig_name));
        }
    }
    
    result.shrink_to_fit();
    
    return result;
}

// private methods

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
    hts_tids_.reserve(hts_header_->n_targets);
    contig_names_.reserve(hts_header_->n_targets);
    
    for (HtsTidType hts_tid {0}; hts_tid < hts_header_->n_targets; ++hts_tid) {
        hts_tids_.emplace(hts_header_->target_name[hts_tid], hts_tid);
        contig_names_.emplace(hts_tid, hts_header_->target_name[hts_tid]);
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

HtslibSamFacade::HtsTidType HtslibSamFacade::get_htslib_tid(const std::string& contig_name) const
{
    return hts_tids_.at(contig_name);
}

const std::string& HtslibSamFacade::contig_name(HtsTidType hts_tid) const
{
    return contig_names_.at(hts_tid);
}

// HtslibIterator
    
auto make_hts_iterator(const hts_idx_t *idx, bam_hdr_t *hdr, const GenomicRegion& region)
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

char extract_base(const uint8_t* hts_sequence, const uint32_t index) noexcept
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
    std::generate_n(std::begin(result), sequence_length,
                    [&i, &hts_sequence] () { return extract_base(hts_sequence, i++); });
    
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
                   [] (const auto op) -> CigarOperation {
                       return CigarOperation {
                           static_cast<CigarOperation::Size>(bam_cigar_oplen(op)),
                           static_cast<CigarOperation::Flag>(bam_cigar_opchr(op))
                       };
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
    result.first_template_segment       = (c.flag & BAM_FREAD1)         != 0;
    result.last_template_segmenet       = (c.flag & BAM_FREAD2)         != 0;
    result.secondary_alignment          = (c.flag & BAM_FSECONDARY)     != 0;
    result.qc_fail                      = (c.flag & BAM_FQCFAIL)        != 0;
    result.duplicate                    = (c.flag & BAM_FDUP)           != 0;
    result.supplementary_alignment      = (c.flag & BAM_FSUPPLEMENTARY) != 0;
    
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

AlignedRead HtslibSamFacade::HtslibIterator::operator*() const
{
    using std::begin; using std::end; using std::next; using std::move;
    
    auto qualities = extract_qualities(hts_bam1_.get());
    
    if (qualities.empty() || qualities[0] == 0xff) {
        throw InvalidBamRecord {hts_facade_.file_path_, extract_read_name(hts_bam1_.get()), "corrupt sequence data"};
    }
    
    auto cigar = extract_cigar_string(hts_bam1_.get());
    const auto& info = hts_bam1_->core;
    auto read_begin_tmp = clipped_begin(cigar, info.pos);
    auto sequence = extract_sequence(hts_bam1_.get());
    
    if (read_begin_tmp < 0) {
        // Then the read hangs off the left of the contig, and we must remove bases, qualities, and
        // adjust the cigar string as we cannot have a negative begin position
        auto overhang_size = std::abs(read_begin_tmp);
        sequence.erase(begin(sequence), next(begin(sequence), overhang_size));
        qualities.erase(begin(qualities), next(begin(qualities), overhang_size));
        
        auto soft_clip_size = cigar.front().size();
        
        if (overhang_size == soft_clip_size) {
            cigar.erase(begin(cigar));
        } else { // then soft_clip_size > overhang_size
            cigar.front() = CigarOperation {soft_clip_size - overhang_size, CigarOperation::Flag::softClipped};
        }
        
        read_begin_tmp = 0;
    }
    
    const auto read_begin = static_cast<AlignedRead::MappingDomain::Position>(read_begin_tmp);
    const auto& contig_name = hts_facade_.contig_name(info.tid);
    
    if (has_multiple_segments(info)) {
        return AlignedRead {
            GenomicRegion {
                contig_name,
                read_begin,
                read_begin + reference_size<AlignedRead::MappingDomain::Position>(cigar)
            },
            move(sequence),
            move(qualities),
            move(cigar),
            mapping_quality(info),
            extract_flags(info),
            hts_facade_.contig_name(info.mtid),
            next_segment_position(info),
            template_length(info),
            extract_next_segment_flags(info)
        };
    } else {
        return AlignedRead {
            GenomicRegion {
                contig_name,
                read_begin,
                read_begin + reference_size<AlignedRead::MappingDomain::Size>(cigar)
            },
            move(sequence),
            move(qualities),
            move(cigar),
            mapping_quality(info),
            extract_flags(info)
        };
    }
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
                       [] (const auto op) {
                           return bam_cigar_oplen(op) > 0;
                       });
}

std::size_t HtslibSamFacade::HtslibIterator::begin() const noexcept
{
    return hts_bam1_->core.pos;
}
    
} // namespace io
} // namespace octopus
