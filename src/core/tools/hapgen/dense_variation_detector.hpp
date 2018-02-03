// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef dense_variation_detector_hpp
#define dense_variation_detector_hpp

#include <vector>

#include "config/common.hpp"
#include "basics/genomic_region.hpp"
#include "core/types/variant.hpp"

namespace octopus { namespace coretools {

class DenseVariationDetector
{
public:
    struct DenseRegion
    {
        enum class RecommendedAction { skip, restrict_lagging };
        GenomicRegion region;
        RecommendedAction action;
    };
    
    DenseVariationDetector() = default;
    
    DenseVariationDetector(double heterozygosity, double heterozygosity_stdev);
    
    DenseVariationDetector(const DenseVariationDetector&)            = default;
    DenseVariationDetector& operator=(const DenseVariationDetector&) = default;
    DenseVariationDetector(DenseVariationDetector&&)                 = default;
    DenseVariationDetector& operator=(DenseVariationDetector&&)      = default;
    
    ~DenseVariationDetector() = default;
    
    std::vector<DenseRegion> detect(const MappableFlatSet<Variant>& variants, const ReadMap& reads) const;

private:
    double expected_heterozygosity_, heterozygosity_stdev_;
};

} // namespace coretools
} // namespace octopus

#endif
