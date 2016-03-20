//
//  downsampler.cpp
//  Octopus
//
//  Created by Daniel Cooke on 19/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "downsampler.hpp"

#include <vector>
#include <deque>
#include <iterator>
#include <algorithm>
#include <random>
#include <cassert>

#include "mappable_algorithms.hpp"
#include "read_utils.hpp"

namespace Octopus
{

bool has_minimum_coverage(const std::vector<unsigned>& required_coverage)
{
    return std::all_of(std::cbegin(required_coverage), std::cend(required_coverage),
                       [] (const auto coverage) { return coverage == 0; });
}

template <typename InputIt>
std::deque<AlignedRead>
sample(const InputIt first_read, const InputIt last_read, const GenomicRegion& region,
       std::vector<unsigned>& required_coverage, const unsigned max_coverage)
{
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::next;
    
    if (first_read == last_read) return {};
    
    const auto positions = decompose(region);
    
    const auto num_positions = positions.size();
    
    assert(num_positions == required_coverage.size());
    
    static std::default_random_engine generator {};
    
    // first pass ensures minimum coverage requirements are satisfied
    
    std::deque<AlignedRead> result {}, unsampled_reads {first_read, last_read};
    
    while (!has_minimum_coverage(required_coverage)) {
        assert(!unsampled_reads.empty());
        
        std::discrete_distribution<unsigned> covers {cbegin(required_coverage), cend(required_coverage)};
        
        const auto sample_position = covers(generator);
        
        const auto overlapped = overlap_range(unsampled_reads, positions[sample_position]);
        
        assert(!overlapped.empty());
        
        const auto num_overlapped = size(overlapped);
        
        std::uniform_int_distribution<std::size_t> read_sampler(0, num_overlapped - 1);
        
        const auto sampled_read_it = std::next(begin(overlapped), read_sampler(generator)).base();
        
        const auto sample_offset        = begin_distance(*sampled_read_it, region);
        const auto num_sample_positions = region_size(*sampled_read_it);
        
        result.emplace_back(std::move(*sampled_read_it));
        unsampled_reads.erase(sampled_read_it);
        
        const auto it = next(begin(required_coverage), sample_offset);
        
        auto it2 = it;
        if (sample_offset + num_sample_positions > num_positions) {
            it2 = end(required_coverage);
        } else {
            it2 = next(it, num_sample_positions);
        }
        
        std::transform(it, it2, it, [] (const auto count) { return (count == 0) ? 0 : count - 1; });
    }
    
    // TODO: second pass increases coverage up to maximum coverage bound
    
    result.shrink_to_fit();
    
    return result;
}

namespace
{
    std::vector<GenomicRegion>
    find_target_coverage_regions(const MappableFlatMultiSet<AlignedRead>& reads,
                                 const unsigned max_coverage, const unsigned min_coverage)
    {
        const auto above_max_coverage_regions = find_high_coverage_regions(reads, max_coverage);
        
        std::vector<GenomicRegion> result {};
        
        if (above_max_coverage_regions.empty()) return result;
        
        result.reserve(above_max_coverage_regions.size());
        
        auto above_min_coverage_regions = find_high_coverage_regions(reads, min_coverage);
        
        std::copy_if(std::make_move_iterator(std::begin(above_min_coverage_regions)),
                     std::make_move_iterator(std::end(above_min_coverage_regions)),
                     std::back_inserter(result),
                     [&above_max_coverage_regions] (const auto& region) {
                         return has_contained(std::cbegin(above_max_coverage_regions),
                                              std::cend(above_max_coverage_regions), region);
                     });
        
        result.shrink_to_fit();
        
        return result;
    }
    
    void add_to_front(std::deque<AlignedRead>& cur_samples, std::deque<AlignedRead>&& new_samples)
    {
        cur_samples.insert(std::begin(cur_samples),
                           std::make_move_iterator(std::begin(new_samples)),
                           std::make_move_iterator(std::end(new_samples)));
    }
    
    template <typename Range>
    auto calculate_minimum_coverages(const Range& reads, const GenomicRegion& region,
                                     const unsigned min_coverage)
    {
        auto result = calculate_positional_coverage(std::cbegin(reads), std::cend(reads), region);
        
        std::transform(std::cbegin(result), std::cend(result), std::begin(result),
                       [min_coverage] (const auto covergae) {
                           return std::min(covergae, min_coverage);
                       });
        
        return result;
    }
} // namespace

std::size_t downsample(MappableFlatMultiSet<AlignedRead>& reads,
                       const unsigned max_coverage, const unsigned min_coverage)
{
    using std::begin; using std::end; using std::make_move_iterator;
    
    const auto sum_reads_before_downsampling = reads.size();
    
    if (sum_reads_before_downsampling == 0) return 0;
    
    const auto regions_to_sample = find_target_coverage_regions(reads, max_coverage, min_coverage);
    
    std::deque<AlignedRead> samples {};
    
    // downsample in reverse order because erasing near back of MappableFlatMultiSet is much
    // cheaper than erasing near front.
    std::for_each(std::make_reverse_iterator(std::cend(regions_to_sample)),
                  std::make_reverse_iterator(std::cbegin(regions_to_sample)),
                  [&] (const auto& region) {
                      //std::cout << "downsampling " << region << std::endl;
                      const auto contained = bases(contained_range(reads, region));
                      
                      auto coverage_requirements = calculate_minimum_coverages(contained, region,
                                                                               min_coverage);
                      
                      add_to_front(samples, sample(begin(contained), end(contained), region,
                                                   coverage_requirements, max_coverage));
                      
                      reads.erase(begin(contained), end(contained));
                  });
    
    reads.insert(make_move_iterator(begin(samples)), make_move_iterator(end(samples)));
    
    return sum_reads_before_downsampling - reads.size();
}

// Downsampler

Downsampler::Downsampler(unsigned max_coverage, unsigned min_coverage)
:
max_coverage_ {max_coverage},
min_coverage_ {min_coverage}
{}

std::size_t Downsampler::downsample(MappableFlatMultiSet<AlignedRead>& reads) const
{
    return Octopus::downsample(reads, max_coverage_, min_coverage_);
}

std::size_t Downsampler::downsample(ReadMap& reads) const
{
    std::size_t result {0};
    
    for (auto& p : reads) {
        result += downsample(p.second);
    }
    
    return result;
}

} // namespace Octopus
