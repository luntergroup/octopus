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
template <typename Region, typename T = unsigned>
class CoverageTracker
{
public:
    using RegionType = Region;
    using DepthType  = T;
    
    CoverageTracker() = default;
    
    CoverageTracker(const CoverageTracker&)            = default;
    CoverageTracker& operator=(const CoverageTracker&) = default;
    CoverageTracker(CoverageTracker&&)                 = default;
    CoverageTracker& operator=(CoverageTracker&&)      = default;
    
    ~CoverageTracker() = default;
    
    template <typename MappableType>
    void add(const MappableType& mappable);
    
    bool any() const noexcept;
    bool any(const Region& region) const noexcept;
    
    std::size_t sum() const noexcept;
    std::size_t sum(const Region& region) const noexcept;
    
    DepthType max() const noexcept;
    DepthType max(const Region& region) const noexcept;
    
    DepthType min() const noexcept;
    DepthType min(const Region& region) const noexcept;
    
    double mean() const noexcept;
    double mean(const Region& region) const noexcept;
    
    double stdev() const noexcept;
    double stdev(const Region& region) const noexcept;
    
    double median() const;
    double median(const Region& region) const;
    
    template <typename OutputIt>
    OutputIt get(const Region& region, OutputIt result) const;
    std::vector<DepthType> get(const Region& region) const;
    
    boost::optional<Region> encompassing_region() const;
    bool is_empty() const noexcept;
    std::size_t num_tracked() const noexcept;
    void clear() noexcept;
    
private:
    std::deque<DepthType> coverage_ = {};
    Region encompassing_region_;
    std::size_t num_tracked_ = 0;
    
    using Iterator     = typename decltype(coverage_)::const_iterator;
    using IteratorPair = std::pair<Iterator, Iterator>;
    
    void do_add(const Region& region);
    IteratorPair range(const Region& region) const;
};

// public methods

template <typename Region, typename T>
template <typename MappableType>
void CoverageTracker<Region, T>::add(const MappableType& mappable)
{
    static_assert(is_region_or_mappable<MappableType>, "MappableType not Mappable");
    do_add(mapped_region(mappable));
}

template <typename Region, typename T>
bool CoverageTracker<Region, T>::any() const noexcept
{
    return std::find_if(std::cbegin(coverage_), std::cend(coverage_),
                        [] (auto depth) noexcept { return depth > 0; }) != std::cend(coverage_);
}

template <typename Region, typename T>
bool CoverageTracker<Region, T>::any(const Region& region) const noexcept
{
    if (octopus::is_empty(region)) return false;
    const auto p = range(region);
    return std::find_if(p.first, p.second, [] (auto depth) noexcept { return depth > 0; }) != p.second;
}

template <typename Region, typename T>
std::size_t CoverageTracker<Region, T>::sum() const noexcept
{
    return std::accumulate(std::cbegin(coverage_), std::cend(coverage_), std::size_t {0});
}

template <typename Region, typename T>
std::size_t CoverageTracker<Region, T>::sum(const Region& region) const noexcept
{
    if (octopus::is_empty(region)) return 0;
    const auto p = range(region);
    return std::accumulate(p.first, p.second, std::size_t {0});
}

template <typename Region, typename T>
T CoverageTracker<Region, T>::max() const noexcept
{
    if (coverage_.empty()) return 0;
    return *std::max_element(std::cbegin(coverage_), std::cend(coverage_));
}

template <typename Region, typename T>
T CoverageTracker<Region, T>::max(const Region& region) const noexcept
{
    if (octopus::is_empty(region)) return 0;
    const auto p = range(region);
    if (p.first == p.second) return 0;
    return *std::max_element(p.first, p.second);
}

template <typename Region, typename T>
T CoverageTracker<Region, T>::min() const noexcept
{
    if (coverage_.empty()) return 0;
    return *std::min_element(std::cbegin(coverage_), std::cend(coverage_));
}

template <typename Region, typename T>
T CoverageTracker<Region, T>::min(const Region& region) const noexcept
{
    if (octopus::is_empty(region)) return 0;
    const auto p = range(region);
    if (p.first == p.second) return 0;
    return *std::min_element(p.first, p.second);
}

template <typename Region, typename T>
double CoverageTracker<Region, T>::mean() const noexcept
{
    if (coverage_.empty()) return 0;
    return maths::mean(coverage_);
}

template <typename Region, typename T>
double CoverageTracker<Region, T>::mean(const Region& region) const noexcept
{
    if (octopus::is_empty(region)) return 0;
    const auto p = range(region);
    return maths::mean(p.first, p.second);
}

template <typename Region, typename T>
double CoverageTracker<Region, T>::stdev() const noexcept
{
    if (coverage_.empty()) return 0;
    return maths::stdev(coverage_);
}

template <typename Region, typename T>
double CoverageTracker<Region, T>::stdev(const Region& region) const noexcept
{
    if (octopus::is_empty(region)) return 0;
    const auto p = range(region);
    return maths::stdev(p.first, p.second);
}

template <typename Region, typename T>
double CoverageTracker<Region, T>::median() const
{
    if (coverage_.empty()) return 0;
    return maths::median(coverage_);
}

template <typename Region, typename T>
double CoverageTracker<Region, T>::median(const Region& region) const
{
    auto range_coverage = this->get(region);
    if (range_coverage.empty()) return 0;
    return maths::median(range_coverage);
}

template <typename Region, typename T>
template <typename OutputIt>
OutputIt CoverageTracker<Region, T>::get(const Region& region, OutputIt result) const
{
    if (coverage_.empty()) {
        return std::fill_n(result, size(region), 0);
    }
    const auto p = range(region);
    if (contains(encompassing_region_, region)) {
        return std::copy(p.first, p.second, result);
    } else {
        using D = typename Region::Distance;
        const auto lhs_pad = std::max(begin_distance(region, encompassing_region_), D {0});
        result = std::fill_n(result, lhs_pad, 0);
        result = std::copy(p.first, p.second, result);
        const auto rhs_pad = std::max(end_distance(encompassing_region_, region), D {0});
        return std::fill_n(result, rhs_pad, 0);
    }
}

template <typename Region, typename T>
std::vector<T> CoverageTracker<Region, T>::get(const Region& region) const
{
    std::vector<T> result(size(region));
    this->get(region, std::begin(result));
    return result;
}

template <typename Region, typename T>
boost::optional<Region> CoverageTracker<Region, T>::encompassing_region() const
{
    if (!is_empty()) {
        return encompassing_region_;
    } else {
        return boost::none;
    }
}

template <typename Region, typename T>
bool CoverageTracker<Region, T>::is_empty() const noexcept
{
    return num_tracked_ == 0;
}

template <typename Region, typename T>
std::size_t CoverageTracker<Region, T>::num_tracked() const noexcept
{
    return num_tracked_;
}

template <typename Region, typename T>
void CoverageTracker<Region, T>::clear() noexcept
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

template <typename Region, typename T>
void CoverageTracker<Region, T>::do_add(const Region& region)
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
        std::transform(first, std::next(first, size(region)), first, [] (auto count) noexcept { return count + 1; });
    }
    ++num_tracked_;
}

template <typename Region, typename T>
typename CoverageTracker<Region, T>::IteratorPair CoverageTracker<Region, T>::range(const Region& region) const
{
    if (coverage_.empty() || !overlaps(region, encompassing_region_)) {
        return {std::end(coverage_), std::end(coverage_)};
    }
    auto range_start_itr = std::begin(coverage_);
    if (begins_before(encompassing_region_, region)) {
        std::advance(range_start_itr, begin_distance(encompassing_region_, region));
    }
    return {range_start_itr, std::next(range_start_itr, overlap_size(region, encompassing_region_))};
}

} // namespace octopus

#endif
