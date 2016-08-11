// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_size_estimator.hpp"

#include <random>
#include <deque>
#include <iterator>
#include <algorithm>

#include "maths.hpp"

namespace octopus {

template <typename ForwardIt, typename RandomGenerator>
ForwardIt random_select(ForwardIt first, ForwardIt last, RandomGenerator& g) {
    std::uniform_int_distribution<std::size_t> dis(0, std::distance(first, last) - 1);
    std::advance(first, dis(g));
    return first;
}

template <typename ForwardIt>
ForwardIt random_select(ForwardIt first, ForwardIt last) {
    static std::default_random_engine gen {};
    return random_select(first, last, gen);
}

auto estimate_read_size(const AlignedRead& read)
{
    return sizeof(AlignedRead)
    // Now the dynamically allocated bits
    + sequence_size(read) * sizeof(char)
    + sequence_size(read) * sizeof(AlignedRead::BaseQuality)
    + read.cigar().size() * sizeof(CigarOperation)
    + contig_name(read).size()
    + (read.has_other_segment() ? sizeof(AlignedRead::Segment) : 0);
}

boost::optional<std::size_t> estimate_mean_read_size(const std::vector<SampleName>& samples,
                                                     const InputRegionMap& input_regions,
                                                     ReadManager& read_manager,
                                                     const unsigned max_sample_size)
{
    assert(!input_regions.empty());
    
    const auto num_samples_per_sample = max_sample_size / samples.size();
    
    std::deque<std::size_t> read_sizes {};
    
    for (const auto& sample : samples) {
        const auto it  = random_select(std::cbegin(input_regions), std::cend(input_regions));
        
        assert(!it->second.empty());
        
        const auto it2 = random_select(std::cbegin(it->second), std::cend(it->second));
        
        auto test_region = read_manager.find_covered_subregion(sample, *it2, num_samples_per_sample);
        
        if (is_empty(test_region)) {
            test_region = expand_rhs(test_region, 1);
        }
        
        const auto reads = read_manager.fetch_reads(sample, test_region);
        
        std::transform(std::cbegin(reads), std::cend(reads), std::back_inserter(read_sizes),
                       estimate_read_size);
    }
    
    if (read_sizes.empty()) {
        return boost::none;
    }
    
    return static_cast<std::size_t>(maths::mean(read_sizes) + maths::stdev(read_sizes));
}

std::size_t default_read_size_estimate() noexcept
{
    return sizeof(AlignedRead) + 300;
}

} // namespace octopus
