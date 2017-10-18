// Copyright (c) 2016 Daniel Cooke
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
#include "core/types/variant.hpp"
#include "io/variant/vcf_record.hpp"
#include "../facets/facet.hpp"
#include "../measures/measure.hpp"

namespace octopus {

class GenomicRegion;
class VcfReader;
class VcfWriter;
class VcfHeader;

namespace csr {

class VariantCallFilter
{
public:
    VariantCallFilter() = delete;
    
    VariantCallFilter(const ReferenceGenome& reference, std::vector<MeasureWrapper> measures);
    
    VariantCallFilter(const VariantCallFilter&)            = delete;
    VariantCallFilter& operator=(const VariantCallFilter&) = delete;
    VariantCallFilter(VariantCallFilter&&)                 = default;
    VariantCallFilter& operator=(VariantCallFilter&&)      = default;
    
    virtual ~VariantCallFilter() = default;
    
    void filter(const VcfReader& source, VcfWriter& dest);
    
protected:
    using MeasureDomain = std::result_of_t<MeasureWrapper(VcfRecord)>;
    using MeasureVector = std::vector<MeasureDomain>;
    
    struct Classification
    {
        enum class Category { filtered, unfiltered } category;
        std::vector<std::string> reasons;
        boost::optional<Phred<double>> quality;
    };
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    std::vector<MeasureWrapper> measures_;
    
private:
    using FacetSet = std::vector<std::string>;
    
    FacetSet facets_;
    
    virtual void annotate(VcfHeader& header) const = 0;
    virtual Classification classify(const MeasureVector& call_measures) const = 0;
    
    VcfRecord filter(const VcfRecord& call) const;
    std::vector<VcfRecord> filter(std::vector<VcfRecord> calls) const;
    FacetSet compute_facets(const VcfRecord& call) const;
    MeasureVector measure(const VcfRecord& call) const;
    void annotate(VcfRecord::Builder& call, Classification status) const;
    void pass(VcfRecord::Builder& call) const;
    void fail(VcfRecord::Builder& call, std::vector<std::string> reasons) const;
};

} // namespace csr

using csr::VariantCallFilter;

} // namespace octopus

#endif
