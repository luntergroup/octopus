// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef variant_call_filter_hpp
#define variant_call_filter_hpp

#include <vector>
#include <string>
#include <cstddef>
#include <type_traits>
#include <functional>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "basics/phred.hpp"
#include "basics/genomic_region.hpp"
#include "core/types/variant.hpp"
#include "io/variant/vcf_header.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_reader.hpp"
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
    };
    
    VariantCallFilter() = delete;
    
    VariantCallFilter(FacetFactory facet_factory,
                      std::vector<MeasureWrapper> measures,
                      OutputOptions output_config);
    
    VariantCallFilter(const VariantCallFilter&)            = delete;
    VariantCallFilter& operator=(const VariantCallFilter&) = delete;
    VariantCallFilter(VariantCallFilter&&)                 = default;
    VariantCallFilter& operator=(VariantCallFilter&&)      = default;
    
    virtual ~VariantCallFilter() = default;
    
    void filter(const VcfReader& source, VcfWriter& dest) const;
    
protected:
    using SampleList    = std::vector<SampleName>;
    using MeasureVector = std::vector<Measure::ResultType>;
    using VcfIterator   = VcfReader::RecordIterator;
    
    struct Classification
    {
        enum class Category { hard_filtered, soft_filtered, unfiltered } category;
        std::vector<std::string> reasons = {};
        boost::optional<Phred<double>> quality = boost::none;
    };
    
    bool can_measure_single_call() const noexcept;
    std::vector<VcfRecord> get_next_block(VcfIterator& first, const VcfIterator& last, const SampleList& samples) const;
    MeasureVector measure(const VcfRecord& call) const;
    std::vector<MeasureVector> measure(const std::vector<VcfRecord>& calls) const;
    void write(const VcfRecord& call, const Classification& classification, VcfWriter& dest) const;
    
private:
    using FacetSet = std::vector<std::string>;
    
    FacetFactory facet_factory_;
    FacetSet facets_;
    std::vector<MeasureWrapper> measures_;
    OutputOptions output_config_;
    
    virtual void annotate(VcfHeader::Builder& header) const = 0;
    virtual void filter(const VcfReader& source, VcfWriter& dest, const SampleList& samples) const = 0;
    
    VcfHeader make_header(const VcfReader& source) const;
    boost::optional<VcfRecord> filter(const VcfRecord& call) const;
    std::vector<VcfRecord> filter(const std::vector<VcfRecord>& calls) const;
    Measure::FacetMap compute_facets(const std::vector<VcfRecord>& calls) const;
    MeasureVector measure(const VcfRecord& call, const Measure::FacetMap& facets) const;
    VcfRecord::Builder construct_template(const VcfRecord& call) const;
    void annotate(VcfRecord::Builder& call, Classification status) const;
    void pass(VcfRecord::Builder& call) const;
    void fail(VcfRecord::Builder& call, std::vector<std::string> reasons) const;
};

} // namespace csr

using csr::VariantCallFilter;

} // namespace octopus

#endif
