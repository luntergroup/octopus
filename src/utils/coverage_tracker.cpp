// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "coverage_tracker.hpp"

#include <iterator>
#include <algorithm>
#include <stdexcept>
#include <cassert>

#include "maths.hpp"

namespace octopus {

// public methods

unsigned CoverageTracker::max_coverage() const noexcept
{
    if (coverage_.empty()) return 0;
    return *std::max_element(std::cbegin(coverage_), std::cend(coverage_));
}
    
unsigned CoverageTracker::max_coverage(const GenomicRegion& region) const noexcept
{
    const auto p = range(region);
    if (p.first == p.second) return 0;
    return *std::max_element(p.first, p.second);
}

unsigned CoverageTracker::min_coverage() const noexcept
{
    if (coverage_.empty()) return 0;
    return *std::min_element(std::cbegin(coverage_), std::cend(coverage_));
}

unsigned CoverageTracker::min_coverage(const GenomicRegion& region) const noexcept
{
    const auto p = range(region);
    if (p.first == p.second) return 0;
    return *std::min_element(p.first, p.second);
}

double CoverageTracker::mean_coverage() const noexcept
{
    if (coverage_.empty()) return 0;
    return maths::mean(coverage_);
}

double CoverageTracker::mean_coverage(const GenomicRegion& region) const noexcept
{
    const auto p = range(region);
    return maths::mean(p.first, p.second);
}

double CoverageTracker::stdev_coverage() const noexcept
{
    if (coverage_.empty()) return 0;
    return maths::stdev(coverage_);
}

double CoverageTracker::stdev_coverage(const GenomicRegion& region) const noexcept
{
    const auto p = range(region);
    return maths::stdev(p.first, p.second);
}

double CoverageTracker::median_coverage(const GenomicRegion& region) const
{
    auto range_coverage = coverage(region);
    if (range_coverage.empty()) return 0;
    const auto first = std::begin(range_coverage);
    const auto nth = std::next(first, range_coverage.size() / 2);
    std::nth_element(first, nth, std::end(range_coverage));
    return *nth;
}

std::vector<unsigned> CoverageTracker::coverage(const GenomicRegion& region) const
{
    const auto p = range(region);
    return std::vector<unsigned> {p.first, p.second};
}
    
GenomicRegion CoverageTracker::encompassing_region() const
{
    return encompassing_region_;
}
    
void CoverageTracker::clear() noexcept
{
    coverage_.clear();
    coverage_.shrink_to_fit();
    num_mappables_added_ = 0;
}

// private methods

void CoverageTracker::do_add(const octopus::GenomicRegion& region)
{
    if (is_empty(region)) return;
    if (num_mappables_added_ == 0) {
        coverage_.assign(size(region), 1);
        encompassing_region_ = region;
    } else {
        if (!is_same_contig(region, encompassing_region_)) {
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
    ++num_mappables_added_;
}

std::pair<CoverageTracker::Iterator, CoverageTracker::Iterator>
CoverageTracker::range(const GenomicRegion& region) const
{
    const auto first = std::next(std::begin(coverage_), begin_distance(encompassing_region_, region));
    return std::make_pair(first, std::next(first, size(region)));
}

} // namespace octopus
