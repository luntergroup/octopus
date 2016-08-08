// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "coverage_tracker.hpp"

#include <iterator>
#include <algorithm>
#include <stdexcept>

#include "maths.hpp"

namespace octopus {

// public methods

unsigned CoverageTracker::max_coverage(const GenomicRegion& region) const
{
    const auto p = range(region);
    if (p.first == p.second) return 0;
    return *std::max_element(p.first, p.second);
}

unsigned CoverageTracker::min_coverage(const GenomicRegion& region) const
{
    const auto p = range(region);
    if (p.first == p.second) return 0;
    return *std::min_element(p.first, p.second);
}

double CoverageTracker::mean_coverage(const GenomicRegion& region) const
{
    const auto p = range(region);
    return maths::mean(p.first, p.second);
}

double CoverageTracker::stdev_coverage(const GenomicRegion& region) const
{
    const auto p = range(region);
    return maths::stdev(p.first, p.second);
}

std::vector<unsigned> CoverageTracker::coverage(const GenomicRegion& region) const
{
    const auto p = range(region);
    return std::vector<unsigned> {p.first, p.second};
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
            encompassing_region_ = encompassing_region(encompassing_region_, region);
        }
        
        const auto first = std::next(std::begin(coverage_), begin_distance(encompassing_region_, region));
        
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
