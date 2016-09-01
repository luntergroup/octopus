// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "variant_call_filter.hpp"

#include <unordered_map>
#include <map>
#include <limits>
#include <numeric>
#include <cmath>

#include <boost/math/distributions/hypergeometric.hpp>

#include "config/common.hpp"
#include "containers/mappable_flat_set.hpp"
#include "containers/mappable_map.hpp"
#include "core/types/variant.hpp"
#include "core/models/haplotype_likelihood_cache.hpp"
#include "basics/genomic_region.hpp"
#include "basics/aligned_read.hpp"
#include "utils/read_stats.hpp"
#include "utils/maths.hpp"
#include "utils/string_utils.hpp"
#include "io/variant/vcf_reader.hpp"
#include "io/variant/vcf_writer.hpp"
#include "io/variant/vcf_header.hpp"
#include "../utils/genotype_reader.hpp"

namespace octopus { namespace csr  {

VariantCallFilter::VariantCallFilter(const ReferenceGenome& reference,
                                     const ReadPipe& read_pipe,
                                     std::vector<MeasureWrapper> measures,
                                     std::size_t max_read_buffer_size)
: measures_ {std::move(measures)}
, reference_ {reference}
, read_pipe_ {read_pipe}
, read_buffer_size_ {max_read_buffer_size}
{}

namespace  {

auto mapped_region(const VcfRecord& call)
{
    const auto begin = call.pos() - 1;
    return GenomicRegion {call.chrom(), begin, begin + static_cast<GenomicRegion::Size>(call.ref().size())};
}

auto mapped_regions(const std::vector<VcfRecord>& calls)
{
    std::vector<GenomicRegion> result {};
    result.reserve(calls.size());
    
    std::transform(std::cbegin(calls), std::cend(calls), std::back_inserter(result),
                   [] (const auto& call) { return mapped_region(call); });
    std::sort(std::begin(result), std::end(result));
    
    return result;
}

auto encompassing_region(const std::vector<VcfRecord>& calls)
{
    return encompassing_region(mapped_regions(calls));
}

auto fetch_reads(const std::vector<VcfRecord>& calls, const ReadPipe& read_pipe)
{
    return read_pipe.fetch_reads(encompassing_region(calls));
}

} // namespace

void VariantCallFilter::register_training_set(const VcfReader& calls, double confidence)
{
    training_sets_.emplace_back(calls, confidence);
}

void VariantCallFilter::filter(const VcfReader& source, VcfWriter& dest)
{
    if (!dest.is_header_written()) {
        auto header = source.fetch_header();
        this->annotate(header);
        dest << header;
    }
    
//    if (this->is_supervised()) {
//        for (const auto& p : training_sets_) {
//            
//        }
//        
//        this->train();
//    }
    
    auto p = source.iterate();
    std::for_each(std::move(p.first), std::move(p.second), [this, &dest] (const VcfRecord& call) {
        auto classification = this->classify(measure(call));
        
        VcfRecord::Builder filtered_call {call};
        
        if (classification.category != Classification::Category::filtered) {
            pass(filtered_call);
        } else {
            fail(filtered_call);
        }
        
        dest.write(filtered_call.build_once());
    });
}

void VariantCallFilter::annotate(VcfRecord::Builder& call) const
{
    
}

VariantCallFilter::MeasureVector VariantCallFilter::measure(const VcfRecord& call) const
{
    return {};
}

void VariantCallFilter::pass(VcfRecord::Builder& call) const
{
    call.set_filter({"FAIL"});
}

void VariantCallFilter::fail(VcfRecord::Builder& call) const
{
    call.set_passed();
}

} // namespace csr
} // namespace octopus
