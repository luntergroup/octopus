// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef dense_variation_detector_hpp
#define dense_variation_detector_hpp

#include <vector>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "basics/genomic_region.hpp"
#include "core/types/variant.hpp"
#include "readpipe/read_pipe.hpp"
#include "utils/input_reads_profiler.hpp"

namespace octopus { namespace coretools {

class DenseVariationDetector
{
public:
    struct Parameters
    {
        enum class Tolerance { low, normal, high };
        double heterozygosity, heterozygosity_stdev;
        Tolerance density_tolerance = Tolerance::normal;
    };
    
    struct DenseRegion
    {
        enum class RecommendedAction { skip, restrict_lagging };
        GenomicRegion region;
        RecommendedAction action;
    };
    
    DenseVariationDetector() = default;
    
    DenseVariationDetector(Parameters params, boost::optional<ReadSetProfile> = boost::none);
    
    DenseVariationDetector(const DenseVariationDetector&)            = default;
    DenseVariationDetector& operator=(const DenseVariationDetector&) = default;
    DenseVariationDetector(DenseVariationDetector&&)                 = default;
    DenseVariationDetector& operator=(DenseVariationDetector&&)      = default;
    
    ~DenseVariationDetector() = default;
    
    std::vector<DenseRegion>
    detect(const MappableFlatSet<Variant>& variants, const ReadMap& reads,
           boost::optional<const ReadPipe::Report&> reads_report = boost::none) const;

private:
    Parameters params_;
    boost::optional<ReadSetProfile> reads_profile_;
    
    double get_max_expected_log_allele_count_per_base() const noexcept;
};

} // namespace coretools
} // namespace octopus

#endif
