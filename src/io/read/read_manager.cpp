// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_manager.hpp"

#include <iterator>
#include <algorithm>
#include <utility>
#include <deque>
#include <numeric>

#include <boost/filesystem/operations.hpp>

#include <basics/aligned_read.hpp>
#include <utils/append.hpp>

#include <iostream> // TEST

namespace octopus {

ReadManager::ReadManager(std::vector<Path> read_file_paths, unsigned max_open_files)
:
max_open_files_ {max_open_files},
num_files_ {static_cast<unsigned>(read_file_paths.size())},
closed_readers_ {std::make_move_iterator(std::begin(read_file_paths)),
                 std::make_move_iterator(std::end(read_file_paths))},
open_readers_ {FileSizeCompare {}},
reader_paths_containing_sample_ {},
possible_regions_in_readers_ {},
samples_ {}
{
    setup_reader_samples_and_regions();
    
    open_initial_files();
    
    samples_.reserve(reader_paths_containing_sample_.size());
    
    for (const auto& pair : reader_paths_containing_sample_) {
        samples_.emplace_back(pair.first);
    }
    
    std::sort(std::begin(samples_), std::end(samples_)); // just for consistency
}

ReadManager::ReadManager(std::initializer_list<Path> read_file_paths)
:
ReadManager {std::vector<Path> {read_file_paths}, static_cast<unsigned>(read_file_paths.size())}
{}

ReadManager::ReadManager(ReadManager&& other)
:
num_files_ {std::move(other.num_files_)}
{
    std::lock_guard<std::mutex> lock {other.mutex_};
    closed_readers_                 = std::move(other.closed_readers_);
    open_readers_                   = std::move(other.open_readers_);
    reader_paths_containing_sample_ = std::move(other.reader_paths_containing_sample_);
    possible_regions_in_readers_    = std::move(other.possible_regions_in_readers_);
    samples_                        = std::move(other.samples_);
}

void swap(ReadManager& lhs, ReadManager& rhs) noexcept
{
    using std::swap;
    if (&lhs == &rhs) return;
    std::lock(lhs.mutex_, rhs.mutex_);
    std::lock_guard<std::mutex> lock_lhs {lhs.mutex_, std::adopt_lock}, lock_rhs {rhs.mutex_, std::adopt_lock};
    //swap(lhs.num_files_, rhs.num_files_);
    swap(lhs.closed_readers_, rhs.closed_readers_);
    swap(lhs.open_readers_, rhs.open_readers_);
    swap(lhs.reader_paths_containing_sample_, rhs.reader_paths_containing_sample_);
    swap(lhs.possible_regions_in_readers_, rhs.possible_regions_in_readers_);
    swap(lhs.samples_, rhs.samples_);
}

bool ReadManager::good() const noexcept
{
    return std::all_of(std::cbegin(open_readers_), std::cend(open_readers_),
                       [] (const auto& p) { return p.second.is_open(); });
}

unsigned ReadManager::num_files() const noexcept
{
    return static_cast<unsigned>(closed_readers_.size() + open_readers_.size());
}

unsigned ReadManager::num_samples() const noexcept
{
    return static_cast<unsigned>(samples_.size());
}

const std::vector<ReadManager::SampleName>& ReadManager::samples() const
{
    return samples_;
}

bool ReadManager::has_reads(const SampleName& sample, const GenomicRegion& region) const
{
    return has_reads({sample}, region);
}

bool ReadManager::has_reads(const std::vector<SampleName>& samples,
                            const GenomicRegion& region) const
{
    auto reader_paths = get_reader_paths_containing_samples(samples);
    
    auto it = partition_open(reader_paths);
    
    while (!reader_paths.empty()) {
        if (std::any_of(it, end(reader_paths),
                        [this, &region, &samples] (const auto& reader_path) {
                            return open_readers_.at(reader_path).has_reads(samples, region);
                        })) {
                            return true;
                        }
        
        reader_paths.erase(it, end(reader_paths));
        it = open_readers(begin(reader_paths), end(reader_paths));
    }
    
    return false;
}

bool ReadManager::has_reads(const GenomicRegion& region) const
{
    auto reader_paths = get_reader_paths_containing_samples(samples());
    
    auto it = partition_open(reader_paths);
    
    while (!reader_paths.empty()) {
        if (std::any_of(it, end(reader_paths),
                        [this, &region] (const auto& reader_path) {
                            return open_readers_.at(reader_path).has_reads(region);
                        })) {
                            return true;
                        }
        
        reader_paths.erase(it, end(reader_paths));
        it = open_readers(begin(reader_paths), end(reader_paths));
    }
    
    return false;
}

std::size_t ReadManager::count_reads(const SampleName& sample, const GenomicRegion& region) const
{
    using std::begin; using std::end; using std::for_each;
    
    std::lock_guard<std::mutex> lock {mutex_};
    
    auto reader_paths = get_possible_reader_paths({sample}, region);
    
    auto it = partition_open(reader_paths);
    
    std::size_t result {0};
    
    while (!reader_paths.empty()) {
        for_each(it, end(reader_paths),
                 [this, &sample, &region, &result] (const auto& reader_path) {
                     result += open_readers_.at(reader_path).count_reads(sample, region);
                 });
        
        reader_paths.erase(it, end(reader_paths));
        it = open_readers(begin(reader_paths), end(reader_paths));
    }
    
    return result;
}

std::size_t ReadManager::count_reads(const std::vector<SampleName>& samples,
                                     const GenomicRegion& region) const
{
    using std::begin; using std::end; using std::for_each;
    
    std::lock_guard<std::mutex> lock {mutex_};
    
    auto reader_paths = get_possible_reader_paths(samples, region);
    
    auto it = partition_open(reader_paths);
    
    std::size_t result {0};
    
    while (!reader_paths.empty()) {
        for_each(it, end(reader_paths),
                 [this, &samples, &region, &result] (const auto& reader_path) {
                     result += open_readers_.at(reader_path).count_reads(samples, region);
                 });
        
        reader_paths.erase(it, end(reader_paths));
        it = open_readers(begin(reader_paths), end(reader_paths));
    }
    
    return result;
}

std::size_t ReadManager::count_reads(const GenomicRegion& region) const
{
    return count_reads(samples(), region);
}

GenomicRegion ReadManager::find_covered_subregion(const SampleName& sample,
                                                  const GenomicRegion& region,
                                                  const std::size_t max_reads) const
{
    return find_covered_subregion(std::vector<SampleName> {sample}, region, max_reads);
}

GenomicRegion ReadManager::find_covered_subregion(const std::vector<SampleName>& samples,
                                                  const GenomicRegion& region,
                                                  const std::size_t max_reads) const
{
    using std::begin; using std::end; using std::next; using std::for_each;
    
    if (samples.empty() || is_empty(region)) return region;
    
    std::lock_guard<std::mutex> lock {mutex_};
    
    auto reader_paths = get_possible_reader_paths(samples, region);
    
    auto it = partition_open(reader_paths);
    
    auto result = head_region(region);
    
    std::deque<unsigned> position_coverage {};
    
    while (!reader_paths.empty()) {
        for_each(it, end(reader_paths),
                 [this, &samples, &region, &result, &position_coverage, max_reads] (const auto& reader_path) {
                     auto p = open_readers_.at(reader_path).find_covered_subregion(samples, region, max_reads);
                     
                     if (is_empty(result) || is_before(p.first, result)) {
                         position_coverage.assign(begin(p.second), end(p.second));
                         result = std::move(p.first);
                         return;
                     }
                     
                     if (is_after(p.first, result)) return;
                     
                     auto overlap_begin = make_pair(begin(position_coverage), begin(p.second));
                     
                     if (begins_before(p.first, result)) {
                         overlap_begin.second = next(begin(p.second), begin_distance(p.first, result));
                         overlap_begin.first = position_coverage.insert(begin(position_coverage),
                                                                        begin(p.second), overlap_begin.second);
                         result = expand_lhs(result, begin_distance(p.first, result));
                     }
                     
                     auto overlap_end = end(position_coverage);
                     
                     if (ends_before(p.first, result)) {
                         overlap_end = std::prev(end(position_coverage), end_distance(p.first, result));
                         overlap_end = position_coverage.erase(overlap_end, end(position_coverage));
                         result = expand_rhs(result, -end_distance(p.first, result));
                     }
                     
                     std::transform(overlap_begin.first, overlap_end, overlap_begin.second,
                                    overlap_begin.first,
                                    [] (const auto curr, const auto x) {
                                        return curr + x;
                                    });
                 });
        
        reader_paths.erase(it, end(reader_paths));
        it = open_readers(begin(reader_paths), end(reader_paths));
    }
    
    if (result == region) return region;
    
    const auto result_begin = std::max(result.begin(), region.begin());
    
    std::partial_sum(begin(position_coverage), end(position_coverage), begin(position_coverage));
    
    const auto limit = std::lower_bound(begin(position_coverage), end(position_coverage), max_reads);
    
    using SizeType = GenomicRegion::Position;
    
    auto result_size = static_cast<SizeType>(std::distance(begin(position_coverage), limit));
    
    if (begins_before(result, region)) {
        result_size -= std::min(result_size, static_cast<SizeType>(begin_distance(result, region)));
    }
    
    const auto result_end = std::min(result_begin + result_size, region.end());
    
    return GenomicRegion {region.contig_name(), result_begin, result_end};
}

GenomicRegion ReadManager::find_covered_subregion(const GenomicRegion& region, const std::size_t max_reads) const
{
    return find_covered_subregion(samples(), region, max_reads);
}

namespace {
    template <typename Container>
    void merge_insert(Container&& src, Container& dst)
    {
        auto it = octopus::utils::append(std::move(src), dst);
        std::inplace_merge(std::begin(dst), it, std::end(dst));
    }
} // namespace

ReadManager::ReadContainer ReadManager::fetch_reads(const SampleName& sample, const GenomicRegion& region) const
{
    using std::begin; using std::end; using std::make_move_iterator; using std::for_each;
    
    std::lock_guard<std::mutex> lock {mutex_};
    
    auto reader_paths = get_possible_reader_paths({sample}, region);
    
    auto it = partition_open(reader_paths);
    
    ReadContainer result {};
    
    while (!reader_paths.empty()) {
        for_each(it, end(reader_paths),
                 [&] (const auto& reader_path) {
                     merge_insert(open_readers_.at(reader_path).fetch_reads(sample, region), result);
                 });
        
        reader_paths.erase(it, end(reader_paths));
        it = open_readers(begin(reader_paths), end(reader_paths));
    }
    
    return result;
}

ReadManager::SampleReadMap ReadManager::fetch_reads(const std::vector<SampleName>& samples,
                                                    const GenomicRegion& region) const
{
    using std::begin; using std::end; using std::make_move_iterator; using std::for_each;
    
    std::lock_guard<std::mutex> lock {mutex_};
    
    auto reader_paths = get_possible_reader_paths(samples, region);
    
    auto it = partition_open(reader_paths);
    
    SampleReadMap result {samples.size()};
    
    for (const auto& sample : samples) {
        result.emplace(std::piecewise_construct, std::forward_as_tuple(sample), std::forward_as_tuple());
    }
    
    while (!reader_paths.empty()) {
        for_each(it, end(reader_paths),
                 [&] (const auto& reader_path) {
                     auto reads = open_readers_.at(reader_path).fetch_reads(samples, region);
                     for (auto&& p : reads) {
                         merge_insert(std::move(p.second), result.at(p.first));
                         p.second.clear();
                         p.second.shrink_to_fit();
                     }
                 });
        
        reader_paths.erase(it, end(reader_paths));
        it = open_readers(begin(reader_paths), end(reader_paths));
    }
    
    return result;
}

ReadManager::SampleReadMap ReadManager::fetch_reads(const GenomicRegion& region) const
{
    return fetch_reads(samples(), region);
}

// Private methods

bool ReadManager::FileSizeCompare::operator()(const Path& lhs, const Path& rhs) const
{
    return boost::filesystem::file_size(lhs) < boost::filesystem::file_size(rhs);
}

void ReadManager::setup_reader_samples_and_regions()
{
    for (const auto& reader_path : closed_readers_) {
        auto reader = make_reader(reader_path);
        add_possible_regions_to_reader_map(reader_path, reader.extract_possible_regions_in_file());
        add_reader_to_sample_map(reader_path, reader.extract_samples());
    }
}

void ReadManager::open_initial_files()
{
    using std::begin; using std::end; using std::cbegin; using std::cend;
    
    std::vector<Path> reader_paths {cbegin(closed_readers_), cend(closed_readers_)};
    
    auto num_files_to_open = std::min(max_open_files_, static_cast<unsigned>(closed_readers_.size()));
    
    std::nth_element(begin(reader_paths), begin(reader_paths) + num_files_to_open, end(reader_paths),
                     FileSizeCompare {});
    
    open_readers(begin(reader_paths), begin(reader_paths) + num_files_to_open);
}

io::ReadReader ReadManager::make_reader(const Path& reader_path) const
{
    return io::ReadReader {reader_path};
}

bool ReadManager::is_open(const Path& reader_path) const noexcept
{
    return open_readers_.count(reader_path) == 1;
}

std::vector<ReadManager::Path>::iterator ReadManager::partition_open(std::vector<Path>& reader_paths) const
{
    return std::partition(std::begin(reader_paths), std::end(reader_paths),
                          [this] (const auto& path) { return !is_open(path); });
}

unsigned ReadManager::num_open_readers() const noexcept
{
    return static_cast<unsigned>(open_readers_.size());
}

unsigned ReadManager::num_reader_spaces() const noexcept
{
    return max_open_files_ - num_open_readers();
}

void ReadManager::open_reader(const Path& reader_path) const
{
    if (num_open_readers() == max_open_files_) { // do we need this?
        close_reader(choose_reader_to_close());
    }
    open_readers_.emplace(reader_path, make_reader(reader_path));
    closed_readers_.erase(reader_path);
}

std::vector<ReadManager::Path>::iterator
ReadManager::open_readers(std::vector<Path>::iterator first, std::vector<Path>::iterator last) const
{
    if (first == last) return first;
    
    auto num_available_spaces = num_reader_spaces();
    auto num_requested_spaces = static_cast<unsigned>(std::distance(first, last));
    
    if (num_requested_spaces <= num_available_spaces) {
        std::for_each(first, last, [this] (const auto& path) { open_reader(path); });
        return first;
    }
    
    auto num_readers_to_close = std::min(num_open_readers(), num_requested_spaces - num_available_spaces);
    
    close_readers(num_readers_to_close);
    
    num_available_spaces += num_readers_to_close;
    
    // partition range so opened readers come last
    auto first_open = std::next(first, num_requested_spaces - num_available_spaces);
    
    std::for_each(first_open, last, [this] (const auto& path) { open_reader(path); });
    
    return first_open;
}

void ReadManager::close_reader(const Path& reader_path) const
{
    // TODO: we can make use of IReadReaderImpl::close to avoid calling deconstructor on the file.
    open_readers_.erase(reader_path);
    closed_readers_.insert(reader_path);
}

ReadManager::Path ReadManager::choose_reader_to_close() const
{
    return open_readers_.begin()->first; // i.e. smallest file size
}

void ReadManager::close_readers(unsigned n) const
{
    for (; n > 0; --n) {
        close_reader(choose_reader_to_close());
    }
}

void ReadManager::add_possible_regions_to_reader_map(const Path& reader_path,
                                                     const std::vector<GenomicRegion>& regions)
{
    for (const auto& region : regions) {
        possible_regions_in_readers_[reader_path][region.contig_name()].emplace(region.contig_region());
    }
}

bool ReadManager::could_reader_contain_region(const Path& reader_path, const GenomicRegion& region) const
{
    if (possible_regions_in_readers_.count(reader_path) == 0) return false;
    if (possible_regions_in_readers_.at(reader_path).count(region.contig_name()) == 0) return false;
    
    return has_overlapped(possible_regions_in_readers_.at(reader_path).at(region.contig_name()),
                          region.contig_region());
}

std::vector<ReadManager::Path>
ReadManager::get_reader_paths_possibly_containing_region(const GenomicRegion& region) const
{
    std::vector<Path> result {};
    result.reserve(num_files_);
    
    for (const auto& reader_path : closed_readers_) {
        if (could_reader_contain_region(reader_path, region)) {
            result.emplace_back(reader_path);
        }
    }
    
    for (const auto& reader : open_readers_) {
        if (could_reader_contain_region(reader.first, region)) {
            result.emplace_back(reader.first);
        }
    }
    
    return result;
}

void ReadManager::add_reader_to_sample_map(const Path& reader_path, const std::vector<SampleName>& samples_in_reader)
{
    for (const auto& sample : samples_in_reader) {
        reader_paths_containing_sample_[sample].emplace_back(reader_path);
    }
}

std::vector<ReadManager::Path>
ReadManager::get_reader_paths_containing_samples(const std::vector<SampleName>& samples) const
{
    std::unordered_set<Path, PathHash> unique_reader_paths {};
    unique_reader_paths.reserve(num_files_);
    
    for (const auto& sample : samples) {
        for (const auto& reader_path : reader_paths_containing_sample_.at(sample)) {
            unique_reader_paths.emplace(reader_path);
        }
    }
    
    return std::vector<Path> {std::begin(unique_reader_paths), std::end(unique_reader_paths)};
}

std::vector<ReadManager::Path>
ReadManager::get_possible_reader_paths(const std::vector<SampleName>& samples, const GenomicRegion& region) const
{
    auto result = get_reader_paths_containing_samples(samples);
    
    auto it = std::remove_if(std::begin(result), std::end(result),
                             [this, &region] (const auto& path) {
                                 return !could_reader_contain_region(path, region);
                             });
    
    result.erase(it, std::end(result));
    
    return result;
}

} // namespace octopus
