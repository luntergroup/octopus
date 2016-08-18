// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef variant_call_filter_hpp
#define variant_call_filter_hpp

#include <vector>
#include <cstddef>
#include <type_traits>

#include <boost/optional.hpp>

#include <config/common.hpp>
#include <basics/phred.hpp>
#include <readpipe/read_pipe.hpp>
#include <core/types/variant.hpp>
#include <io/reference/reference_genome.hpp>
#include <io/variant/vcf_record.hpp>

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
    
    VariantCallFilter(const ReferenceGenome& reference,
                      const ReadPipe& read_pipe,
                      std::vector<MeasureWrapper> measures,
                      std::size_t max_read_buffer_size);
    
    VariantCallFilter(const VariantCallFilter&)            = delete;
    VariantCallFilter& operator=(const VariantCallFilter&) = delete;
    VariantCallFilter(VariantCallFilter&&)                 = default;
    VariantCallFilter& operator=(VariantCallFilter&&)      = default;
    
    virtual ~VariantCallFilter() = default;
    
    virtual bool is_supervised() const noexcept { return false; }
    
    void register_training_set(const VcfReader& calls, double confidence);
    
    void filter(const VcfReader& source, VcfWriter& dest);
    
protected:
    using MeasureDomain = std::result_of_t<MeasureWrapper(VcfRecord)>;
    using MeasureVector = std::vector<MeasureDomain>;
    
    struct Classification
    {
        enum class Category { filtered, unfiltered } category;
        boost::optional<Phred<double>> quality;
    };
    
    std::vector<MeasureWrapper> measures_;
    
private:
    const ReferenceGenome& reference_;
    const ReadPipe& read_pipe_;
    
    std::size_t read_buffer_size_;
    
    std::deque<std::pair<std::reference_wrapper<const VcfReader>, double>> training_sets_;
    
    virtual void annotate(VcfHeader& header) const = 0;
    virtual void register_training_point(const MeasureVector& call_measures, double confidence) {};
    virtual void train() {};
    virtual Classification classify(const MeasureVector& call_measures) const = 0;
    
    void annotate(VcfRecord::Builder& call) const;
    MeasureVector measure(const VcfRecord& call) const;
    void pass(VcfRecord::Builder& call) const;
    void fail(VcfRecord::Builder& call) const;
};

class VariantCallFilterWrapper
{
public:
    VariantCallFilterWrapper() = delete;
    
    VariantCallFilterWrapper(std::unique_ptr<VariantCallFilter> filter);
    
    ~VariantCallFilterWrapper() = default;
    
    bool is_supervised() const noexcept;
    
    void register_truth(const VcfReader& calls);
    
    void filter(const VcfReader& source, VcfWriter& dest) const;
private:
    std::unique_ptr<VariantCallFilter> filter_;
};

} // namespace csr
} // namespace octopus

#endif /* variant_call_filter_hpp */
