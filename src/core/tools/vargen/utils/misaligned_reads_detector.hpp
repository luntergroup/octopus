// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef misaligned_reads_detector_hpp
#define misaligned_reads_detector_hpp

#include <vector>
#include <unordered_map>
#include <functional>
#include <cmath>

#include "config/common.hpp"
#include "basics/genomic_region.hpp"
#include "basics/aligned_read.hpp"
#include "utils/coverage_tracker.hpp"
#include "io/reference/reference_genome.hpp"

namespace octopus { namespace coretools {

class MisalignedReadsDetector
{
public:
    struct Options
    {
        AlignedRead::BaseQuality snv_threshold;
        double snv_penalty = 1, indel_penalty = 1, clip_penalty = 1;
        double max_expected_mutation_rate = 1e-3;
        double min_ln_prob_correctly_aligned = std::log(0.0001);
        unsigned max_unpenalised_clip_size = 3;
    };
    
    MisalignedReadsDetector() = delete;
    
    MisalignedReadsDetector(const ReferenceGenome& reference);
    MisalignedReadsDetector(const ReferenceGenome& reference, Options options);
    
    MisalignedReadsDetector(const MisalignedReadsDetector&)            = default;
    MisalignedReadsDetector& operator=(const MisalignedReadsDetector&) = default;
    MisalignedReadsDetector(MisalignedReadsDetector&&)                 = default;
    MisalignedReadsDetector& operator=(MisalignedReadsDetector&&)      = default;
    
    ~MisalignedReadsDetector() = default;
    
    void add(const SampleName& sample, const AlignedRead& read);
    template <typename ForwardIterator>
    void add(const SampleName& sample, ForwardIterator first_read, ForwardIterator last_read);
    
    std::vector<GenomicRegion> generate(const GenomicRegion& region) const;
    
    void clear() noexcept;
    
private:
    using CoverageTrackerMap = std::unordered_map<SampleName, CoverageTracker<GenomicRegion>>;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    Options options_;
    CoverageTrackerMap coverage_tracker_, likely_misaligned_coverage_tracker_;
    
    bool is_likely_misaligned(const AlignedRead& read) const;
};

template <typename ForwardIterator>
void MisalignedReadsDetector::add(const SampleName& sample, ForwardIterator first_read, ForwardIterator last_read)
{
    auto& coverage_tracker = coverage_tracker_[sample];
    auto& likely_misaligned_coverage_tracker = likely_misaligned_coverage_tracker_[sample];
    std::for_each(first_read, last_read, [&] (const AlignedRead& read) {
        coverage_tracker.add(read);
        if (is_likely_misaligned(read)) {
            likely_misaligned_coverage_tracker.add(read);
        }
    });
}

} // namespace coretools
} // namespace octopus

#endif
