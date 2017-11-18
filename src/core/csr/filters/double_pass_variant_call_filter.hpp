// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef double_variant_call_filter_hpp
#define double_variant_call_filter_hpp

#include <vector>

#include <boost/optional.hpp>

#include "logging/progress_meter.hpp"
#include "basics/genomic_region.hpp"
#include "logging/logging.hpp"
#include "variant_call_filter.hpp"

namespace octopus { namespace csr {

class DoublePassVariantCallFilter : public VariantCallFilter
{
public:
    DoublePassVariantCallFilter() = delete;
    
    DoublePassVariantCallFilter(FacetFactory facet_factory,
                                std::vector<MeasureWrapper> measures,
                                OutputOptions output_config,
                                boost::optional<ProgressMeter&> progress);
    
    DoublePassVariantCallFilter(const DoublePassVariantCallFilter&)            = delete;
    DoublePassVariantCallFilter& operator=(const DoublePassVariantCallFilter&) = delete;
    DoublePassVariantCallFilter(DoublePassVariantCallFilter&&)                 = default;
    DoublePassVariantCallFilter& operator=(DoublePassVariantCallFilter&&)      = default;
    
    virtual ~DoublePassVariantCallFilter() override = default;
    
private:
    mutable boost::optional<logging::InfoLogger> info_log_;
    mutable boost::optional<ProgressMeter&> progress_;
    mutable boost::optional<GenomicRegion::ContigName> current_contig_;
    
    virtual void record(std::size_t call_idx, const MeasureVector& measures) const = 0;
    virtual Classification classify(std::size_t call_idx) const = 0;
    
    void filter(const VcfReader& source, VcfWriter& dest, const SampleList& samples) const override;
    
    void make_registration_pass(const VcfReader& source, const SampleList& samples) const;
    void record(const VcfRecord& call, std::size_t idx) const;
    void record(const std::vector<VcfRecord>& calls, std::size_t first_idx) const;
    void make_filter_pass(const VcfReader& source, VcfWriter& dest) const;
    void filter(const VcfRecord& call, std::size_t idx, VcfWriter& dest) const;
    void log_progress(const GenomicRegion& region) const;
};

} // namespace csr
} // namespace octopus

#endif
