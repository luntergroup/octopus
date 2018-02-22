// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef coverage_tracker_hpp
#define coverage_tracker_hpp

#include <vector>
#include <deque>
#include <cstddef>
#include <utility>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cassert>

#include <boost/optional.hpp>

#include "concepts/mappable.hpp"
#include "maths.hpp"

namespace octopus {

/**
 CoverageTracker provides an efficient method for tracking coverage statistics over a range
 of Mappable objects without having to store the entire collection.
 */
template <typename Region>
class CoverageTracker
{
public:
    CoverageTracker() = default;
    
    CoverageTracker(const CoverageTracker&)            = default;
    CoverageTracker& operator=(const CoverageTracker&) = default;
    CoverageTracker(CoverageTracker&&)                 = default;
    CoverageTracker& operator=(CoverageTracker&&)      = default;
    
    ~CoverageTracker() = default;
    
    template <typename MappableType>
    void add(const MappableType& mappable);
    
    std::size_t total_coverage() const noexcept;
    std::size_t total_coverage(const Region& region) const noexcept;
    
    unsigned max_coverage() const noexcept;
    unsigned max_coverage(const Region& region) const noexcept;
    
    unsigned min_coverage() const noexcept;
    unsigned min_coverage(const Region& region) const noexcept;
    
    double mean_coverage() const noexcept;
    double mean_coverage(const Region& region) const noexcept;
    
    double stdev_coverage() const noexcept;
    double stdev_coverage(const Region& region) const noexcept;
    
    double median_coverage(const Region& region) const;
    
    std::vector<unsigned> coverage(const Region& region) const;
    
    boost::optional<Region> encompassing_region() const;
    
    bool is_empty() const noexcept;
    
    std::size_t num_tracked() const noexcept;
    
    void clear() noexcept;
    
private:
    std::deque<unsigned> coverage_ = {};
    Region encompassing_region_;
    std::size_t num_tracked_ = 0;
    
    using Iterator = typename decltype(coverage_)::const_iterator;
    
    void do_add(const Region& region);
    std::pair<Iterator, Iterator> range(const Region& region) const;
};

// non-member methods

template <typename Region>
std::vector<Region> get_covered_regions(const CoverageTracker<Region>& tracker, const Region& region)
{
    const auto depths = tracker.coverage(region);
    return select_regions(region, depths, [] (unsigned depth) { return depth > 0; });
}

template <typename Region>
std::vector<Region> get_covered_regions(const CoverageTracker<Region>& tracker)
{
    const auto tracker_region = tracker.encompassing_region();
    if (tracker_region) {
        return get_covered_regions(tracker, *tracker_region);
    } else {
        return {};
    }
}

// public methods

template <typename Region>
template <typename MappableType>
void CoverageTracker<Region>::add(const MappableType& mappable)
{
    static_assert(is_region_or_mappable<MappableType>, "MappableType not Mappable");
    do_add(mapped_region(mappable));
}

template <typename Region>
std::size_t CoverageTracker<Region>::total_coverage() const noexcept
{
    return std::accumulate(std::cbegin(coverage_), std::cend(coverage_), std::size_t {0});
}

template <typename Region>
std::size_t CoverageTracker<Region>::total_coverage(const Region& region) const noexcept
{
    if (octopus::is_empty(region)) return 0;
    const auto p = range(region);
    if (p.first == p.second) return 0;
    return std::accumulate(p.first, p.second, std::size_t {0});
}

template <typename Region>
unsigned CoverageTracker<Region>::max_coverage() const noexcept
{
    if (coverage_.empty()) return 0;
    return *std::max_element(std::cbegin(coverage_), std::cend(coverage_));
}

template <typename Region>
unsigned CoverageTracker<Region>::max_coverage(const Region& region) const noexcept
{
    if (octopus::is_empty(region)) return 0;
    const auto p = range(region);
    if (p.first == p.second) return 0;
    return *std::max_element(p.first, p.second);
}

template <typename Region>
unsigned CoverageTracker<Region>::min_coverage() const noexcept
{
    if (coverage_.empty()) return 0;
    return *std::min_element(std::cbegin(coverage_), std::cend(coverage_));
}

template <typename Region>
unsigned CoverageTracker<Region>::min_coverage(const Region& region) const noexcept
{
    if (octopus::is_empty(region)) return 0;
    const auto p = range(region);
    if (p.first == p.second) return 0;
    return *std::min_element(p.first, p.second);
}

template <typename Region>
double CoverageTracker<Region>::mean_coverage() const noexcept
{
    if (coverage_.empty()) return 0;
    return maths::mean(coverage_);
}

template <typename Region>
double CoverageTracker<Region>::mean_coverage(const Region& region) const noexcept
{
    if (octopus::is_empty(region)) return 0;
    const auto p = range(region);
    return maths::mean(p.first, p.second);
}

template <typename Region>
double CoverageTracker<Region>::stdev_coverage() const noexcept
{
    if (coverage_.empty()) return 0;
    return maths::stdev(coverage_);
}

template <typename Region>
double CoverageTracker<Region>::stdev_coverage(const Region& region) const noexcept
{
    if (octopus::is_empty(region)) return 0;
    const auto p = range(region);
    return maths::stdev(p.first, p.second);
}

template <typename Region>
double CoverageTracker<Region>::median_coverage(const Region& region) const
{
    auto range_coverage = coverage(region);
    if (range_coverage.empty()) return 0;
    const auto first = std::begin(range_coverage);
    const auto nth = std::next(first, range_coverage.size() / 2);
    std::nth_element(first, nth, std::end(range_coverage));
    return *nth;
}

template <typename Region>
std::vector<unsigned> CoverageTracker<Region>::coverage(const Region& region) const
{
    if (coverage_.empty()) return std::vector<unsigned>(size(region), 0);
    const auto p = range(region);
    if (!contains(encompassing_region_, region)) {
        std::vector<unsigned> result(size(region), 0);
        const auto d = std::max(begin_distance(region, encompassing_region_), GenomicRegion::Distance {0});
        std::copy(p.first, p.second, std::next(std::begin(result), d));
        return result;
    }
    return std::vector<unsigned> {p.first, p.second};
}

template <typename Region>
boost::optional<Region> CoverageTracker<Region>::encompassing_region() const
{
    if (!is_empty()) {
        return encompassing_region_;
    } else {
        return boost::none;
    }
}

template <typename Region>
bool CoverageTracker<Region>::is_empty() const noexcept
{
    return num_tracked_ == 0;
}

template <typename Region>
std::size_t CoverageTracker<Region>::num_tracked() const noexcept
{
    return num_tracked_;
}

template <typename Region>
void CoverageTracker<Region>::clear() noexcept
{
    coverage_.clear();
    coverage_.shrink_to_fit();
    num_tracked_ = 0;
}

// private methods

namespace detail {

inline bool is_same_contig_helper(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return true;
}

inline bool is_same_contig_helper(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return is_same_contig(lhs, rhs);
}

} // namespace detail

template <typename Region>
void CoverageTracker<Region>::do_add(const Region& region)
{
    if (octopus::is_empty(region)) return;
    if (num_tracked_ == 0) {
        coverage_.assign(size(region), 1);
        encompassing_region_ = region;
    } else {
        if (!detail::is_same_contig_helper(region, encompassing_region_)) {
            throw std::runtime_error {"CoverageTracker: contig mismatch"};
        }
        bool region_change {false};
        if (begins_before(region, encompassing_region_)) {
            coverage_.insert(std::cbegin(coverage_), left_overhang_size(region, encompassing_region_), 0);
            region_change = true;
        }
        if (ends_before(encompassing_region_, region)) {
            coverage_.insert(std::cend(coverage_), right_overhang_size(region, encompassing_region_), 0);
            region_change = true;
        }
        if (region_change) {
            encompassing_region_ = octopus::encompassing_region(encompassing_region_, region);
        }
        const auto first = std::next(std::begin(coverage_), begin_distance(encompassing_region_, region));
        assert(first < std::end(coverage_));
        assert(std::next(first, size(region)) <= std::end(coverage_));
        std::transform(first, std::next(first, size(region)), first, [] (auto count) { return count + 1; });
    }
    ++num_tracked_;
}

template <typename Region>
std::pair<typename CoverageTracker<Region>::Iterator, typename CoverageTracker<Region>::Iterator>
CoverageTracker<Region>::range(const Region& region) const
{
    if (coverage_.empty() || !overlaps(region, encompassing_region_)) {
        return {std::end(coverage_), std::end(coverage_)};
    }
    auto first = std::begin(coverage_);
    if (begins_before(encompassing_region_, region)) {
        std::advance(first, begin_distance(encompassing_region_, region));
    }
    return {first, std::next(first, overlap_size(region, encompassing_region_))};
}

} // namespace octopus

#endif
