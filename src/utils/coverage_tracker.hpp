// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef coverage_tracker_hpp
#define coverage_tracker_hpp

#include <vector>
#include <deque>
#include <cstddef>
#include <utility>

#include <concepts/mappable.hpp>
#include <basics/genomic_region.hpp>

namespace octopus {

/**
 CoverageTracker provides an efficient method for tracking coverage statistics over a range
 of mappable objects without having to store the entire collection.
 */
class CoverageTracker
{
public:
    CoverageTracker() = default;
    
    CoverageTracker(const CoverageTracker&)            = default;
    CoverageTracker& operator=(const CoverageTracker&) = default;
    CoverageTracker(CoverageTracker&&)                 = default;
    CoverageTracker& operator=(CoverageTracker&&)      = default;
    
    ~CoverageTracker() = default;
    
    template <typename T>
    void add(const T& mappable);
    
    unsigned max_coverage(const GenomicRegion& region) const;
    
    unsigned min_coverage(const GenomicRegion& region) const;
    
    double mean_coverage(const GenomicRegion& region) const;
    
    double stdev_coverage(const GenomicRegion& region) const;
    
    std::vector<unsigned> coverage(const GenomicRegion& region) const;
    
    void clear() noexcept;
    
private:
    std::deque<unsigned> coverage_;
    
    GenomicRegion encompassing_region_;
    
    std::size_t num_mappables_added_;
    
    void do_add(const GenomicRegion& region);
    
    using Iterator = decltype(coverage_)::const_iterator;
    
    std::pair<Iterator, Iterator> range(const GenomicRegion& region) const;
};

template <typename T>
void CoverageTracker::add(const T& mappable)
{
    static_assert(is_region_or_mappable<T>, "T must be a region or mappable");
    do_add(mapped_region(mappable));
}

} // namespace octopus

#endif
