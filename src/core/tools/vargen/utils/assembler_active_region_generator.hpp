// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef assembler_active_region_generator_hpp
#define assembler_active_region_generator_hpp

#include <vector>
#include <unordered_map>
#include <functional>

#include "config/common.hpp"
#include "basics/aligned_read.hpp"
#include "utils/coverage_tracker.hpp"
#include "io/reference/reference_genome.hpp"

namespace octopus {

class GenomicRegion;

namespace coretools {

class AssemblerActiveRegionGenerator
{
public:
    struct Options
    {
        enum class TriggerType { snv, indel, structual };
        // TriggerType::snv looks for reads containing SNVs
        // TriggerType::indel looks for reads containing indels and soft clipping
        // TriggerType::structual looks for deletion hotspots
        std::vector<TriggerType> trigger_types = {TriggerType::indel};
        AlignedRead::BaseQuality trigger_quality = 10;
        AlignedRead::MappingDomain::Size trigger_clip_size = 2;
        double min_expected_mutation_frequency = 0.1;
    };
    
    AssemblerActiveRegionGenerator() = delete;
    
    AssemblerActiveRegionGenerator(const ReferenceGenome& reference);
    AssemblerActiveRegionGenerator(const ReferenceGenome& reference, Options options);
    
    AssemblerActiveRegionGenerator(const AssemblerActiveRegionGenerator&)            = default;
    AssemblerActiveRegionGenerator& operator=(const AssemblerActiveRegionGenerator&) = default;
    AssemblerActiveRegionGenerator(AssemblerActiveRegionGenerator&&)                 = default;
    AssemblerActiveRegionGenerator& operator=(AssemblerActiveRegionGenerator&&)      = default;
    
    ~AssemblerActiveRegionGenerator() = default;
    
    void add(const SampleName& sample, const AlignedRead& read);
    
    template <typename ForwardIterator>
    void add(const SampleName& sample, ForwardIterator first_read, ForwardIterator last_read);
    
    std::vector<GenomicRegion> generate(const GenomicRegion& region) const;

    void clear() noexcept;
    
private:
    using CoverageTrackerMap = std::unordered_map<SampleName, CoverageTracker<GenomicRegion>>;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    bool snvs_interesting_ = false, indels_interesting_ = true, structual_interesting_ = false;
    AlignedRead::BaseQuality trigger_quality_ = 10;
    AlignedRead::MappingDomain::Size trigger_clip_size_ = 2;
    double min_expected_mutation_frequency_;
    CoverageTrackerMap coverage_tracker_, interesting_read_coverages_, clipped_coverage_tracker_;
    
    bool is_interesting(const AlignedRead& read) const;
};

template <typename ForwardIterator>
void AssemblerActiveRegionGenerator::add(const SampleName& sample, ForwardIterator first_read, ForwardIterator last_read)
{
    auto& coverage_tracker = coverage_tracker_[sample];
    auto& interesting_coverage_tracker = interesting_read_coverages_[sample];
    std::for_each(first_read, last_read, [&] (const AlignedRead& read) {
        coverage_tracker.add(read);
        if (is_interesting(read)) {
            interesting_coverage_tracker.add(read);
        }
    });
    if (structual_interesting_) {
        auto& clipped_coverage_tracker = clipped_coverage_tracker_[sample];
        std::for_each(first_read, last_read, [&] (const AlignedRead& read) {
            if (is_soft_clipped(read)) {
                clipped_coverage_tracker.add(clipped_mapped_region(read));
            } else {
                clipped_coverage_tracker.add(read);
            }
        });
    }
}
    
} // namespace coretools
} // namespace octopus

#endif

