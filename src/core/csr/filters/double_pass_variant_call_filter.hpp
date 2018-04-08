// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef double_variant_call_filter_hpp
#define double_variant_call_filter_hpp

#include <vector>
#include <cstddef>

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
                                ConcurrencyPolicy threading,
                                boost::optional<ProgressMeter&> progress);
    
    DoublePassVariantCallFilter(const DoublePassVariantCallFilter&)            = delete;
    DoublePassVariantCallFilter& operator=(const DoublePassVariantCallFilter&) = delete;
    DoublePassVariantCallFilter(DoublePassVariantCallFilter&&)                 = default;
    DoublePassVariantCallFilter& operator=(DoublePassVariantCallFilter&&)      = default;
    
    virtual ~DoublePassVariantCallFilter() override = default;

protected:
    using Log = logging::InfoLogger;
    
private:
    mutable boost::optional<Log> info_log_;
    mutable boost::optional<ProgressMeter&> progress_;
    mutable boost::optional<GenomicRegion::ContigName> current_contig_;
    
    virtual void log_registration_pass(Log& log) const;
    virtual void prepare_for_registration(const SampleList& samples) const {};
    virtual void record(std::size_t call_idx, std::size_t sample_idx, MeasureVector measures) const = 0;
    virtual void prepare_for_classification(boost::optional<Log>& log) const = 0;
    virtual void log_filter_pass_start(Log& log) const;
    virtual Classification classify(std::size_t call_idx, std::size_t sample_idx) const = 0;
    virtual Classification merge(const std::vector<Classification>& sample_classifications) const;
    
    void filter(const VcfReader& source, VcfWriter& dest, const SampleList& samples) const override;
    
    void make_registration_pass(const VcfReader& source, const SampleList& samples) const;
    void record(const VcfRecord& call, std::size_t record_idx, const SampleList& samples) const;
    void record(const CallBlock& block, std::size_t record_idx, const SampleList& samples) const;
    void record(const std::vector<CallBlock>& blocks, std::size_t record_idx, const SampleList& samples) const;
    void record(const VcfRecord& call, const MeasureVector& measures, std::size_t record_idx, const SampleList& samples) const;
    void record(const CallBlock& block, const MeasureBlock& measures, std::size_t record_idx, const SampleList& samples) const;
    void make_filter_pass(const VcfReader& source, const SampleList& samples, VcfWriter& dest) const;
    std::vector<Classification> classify(std::size_t call_idx, const SampleList& samples) const;
    void filter(const VcfRecord& call, std::size_t idx, const SampleList& samples, VcfWriter& dest) const;
    void log_progress(const GenomicRegion& region) const;
};

} // namespace csr
} // namespace octopus

#endif
