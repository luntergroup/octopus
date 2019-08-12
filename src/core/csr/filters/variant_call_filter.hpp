// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef variant_call_filter_hpp
#define variant_call_filter_hpp

#include <vector>
#include <string>
#include <cstddef>
#include <type_traits>
#include <functional>
#include <future>
#include <unordered_set>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "basics/phred.hpp"
#include "basics/genomic_region.hpp"
#include "core/types/variant.hpp"
#include "io/variant/vcf_header.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_reader.hpp"
#include "utils/thread_pool.hpp"
#include "logging/logging.hpp"
#include "../facets/facet.hpp"
#include "../facets/facet_factory.hpp"
#include "../measures/measure.hpp"

namespace octopus {

class VcfWriter;
class VcfHeader;

namespace csr {

class VariantCallFilter
{
public:
    struct OutputOptions
    {
        bool emit_sites_only = false;
        bool clear_existing_filters = true;
        bool clear_info = false;
        bool annotate_all_active_measures = false;
        std::unordered_set<std::string> annotations = {};
    };
    
    struct ConcurrencyPolicy
    {
        boost::optional<unsigned> max_threads = boost::none;
    };
    
    VariantCallFilter() = delete;
    
    VariantCallFilter(FacetFactory facet_factory,
                      std::vector<MeasureWrapper> measures,
                      OutputOptions output_config,
                      ConcurrencyPolicy threading);
    
    VariantCallFilter(const VariantCallFilter&)            = delete;
    VariantCallFilter& operator=(const VariantCallFilter&) = delete;
    VariantCallFilter(VariantCallFilter&&)                 = delete;
    VariantCallFilter& operator=(VariantCallFilter&&)      = delete;
    
    virtual ~VariantCallFilter() = default;
    
    std::string name() const;
    
    void filter(const VcfReader& source, VcfWriter& dest) const;
    
protected:
    using SampleList    = std::vector<SampleName>;
    using MeasureVector = std::vector<Measure::ResultType>;
    using VcfIterator   = VcfReader::RecordIterator;
    using CallBlock     = std::vector<VcfRecord>;
    using MeasureBlock  = std::vector<MeasureVector>;
    
    struct Classification
    {
        enum class Category { hard_filtered, soft_filtered, unfiltered } category;
        std::vector<std::string> reasons = {};
        boost::optional<Phred<double>> quality = boost::none;
    };
    using ClassificationList = std::vector<Classification>;
    
    std::vector<MeasureWrapper> measures_;
    mutable boost::optional<logging::DebugLogger> debug_log_;
    
    virtual Classification merge(const ClassificationList& sample_classifications, const MeasureVector& measures) const;
    virtual Classification merge(const ClassificationList& sample_classifications) const;
    
    bool can_measure_single_call() const noexcept;
    bool can_measure_multiple_blocks() const noexcept;
    CallBlock read_next_block(VcfIterator& first, const VcfIterator& last, const SampleList& samples) const;
    std::vector<CallBlock> read_next_blocks(VcfIterator& first, const VcfIterator& last, const SampleList& samples) const;
    MeasureVector measure(const VcfRecord& call) const;
    MeasureBlock measure(const CallBlock& block) const;
    std::vector<MeasureBlock> measure(const std::vector<CallBlock>& blocks) const;
    void write(const VcfRecord& call, const Classification& classification, VcfWriter& dest) const;
    void write(const VcfRecord& call, const Classification& classification,
               const SampleList& samples, const ClassificationList& sample_classifications,
               VcfWriter& dest) const;
    bool measure_annotations_requested() const noexcept;
    void annotate(VcfRecord::Builder& call, const MeasureVector& measures, const VcfHeader& header) const;
    Phred<double> compute_joint_probability(const std::vector<Phred<double>>& qualities) const;
    std::vector<std::string> compute_reason_union(const ClassificationList& sample_classifications) const;
    
private:
    using FacetNameSet = std::vector<std::string>;
    
    FacetFactory facet_factory_;
    FacetNameSet facet_names_;
    OutputOptions output_config_;
    std::vector<MeasureWrapper> duplicate_measures_;
    
    mutable ThreadPool workers_;
    
    virtual std::string do_name() const = 0;
    virtual void annotate(VcfHeader::Builder& header) const = 0;
    virtual void filter(const VcfReader& source, VcfWriter& dest, const VcfHeader& dest_header) const = 0;
    virtual boost::optional<std::string> call_quality_name() const { return boost::none; }
    virtual boost::optional<std::string> genotype_quality_name() const { return boost::none; }
    virtual boost::optional<Phred<double>> compute_joint_quality(const ClassificationList& sample_classifications, const MeasureVector& measures) const;
    virtual bool is_soft_filtered(const ClassificationList& sample_classifications, boost::optional<Phred<double>> joint_quality,
                                  const MeasureVector& measures, std::vector<std::string>& reasons) const;
    
    VcfHeader make_header(const VcfReader& source) const;
    Measure::FacetMap compute_facets(const CallBlock& block) const;
    std::vector<Measure::FacetMap> compute_facets(const std::vector<CallBlock>& blocks) const;
    MeasureBlock measure(const CallBlock& block, const Measure::FacetMap& facets) const;
    MeasureVector measure(const VcfRecord& call, const Measure::FacetMap& facets) const;
    VcfRecord::Builder construct_template(const VcfRecord& call) const;
    bool is_requested_annotation(const MeasureWrapper& measure) const noexcept;
    bool is_hard_filtered(const Classification& classification) const noexcept;
    void annotate(VcfRecord::Builder& call, const SampleList& samples, const ClassificationList& sample_classifications) const;
    void annotate(VcfRecord::Builder& call, const SampleName& sample, Classification status) const;
    void annotate(VcfRecord::Builder& call, Classification status) const;
    void pass(const SampleName& sample, VcfRecord::Builder& call) const;
    void pass(VcfRecord::Builder& call) const;
    void fail(const SampleName& sample, VcfRecord::Builder& call, std::vector<std::string> reasons) const;
    void fail(VcfRecord::Builder& call, std::vector<std::string> reasons) const;
    bool is_multithreaded() const noexcept;
    unsigned max_concurrent_blocks() const noexcept;
};

} // namespace csr

using csr::VariantCallFilter;

} // namespace octopus

#endif
