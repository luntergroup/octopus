//
//  variant_call_filter.cpp
//  Octopus
//
//  Created by Daniel Cooke on 31/05/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "variant_call_filter.hpp"

#include "vcf_header.hpp"
#include "vcf_record.hpp"

#include "aligned_read.hpp"
#include "read_utils.hpp"

#include "haplotype_likelihood_cache.hpp"

namespace Octopus
{
    VariantCallFilter::VariantCallFilter(ReadManager& read_manager)
    :
    read_manager_ {read_manager}
    {}
    
    namespace
    {
        auto mapped_region(const VcfRecord& call)
        {
            using SizeType = GenomicRegion::SizeType;
            const auto begin = call.pos() - 1;
            return GenomicRegion {call.chrom(), begin, begin + static_cast<SizeType>(call.ref().size())};
        }
    } // namespace
    
    void VariantCallFilter::filter(VcfReader& source, VcfWriter& dest)
    {
        if (!dest.is_header_written()) {
            dest.write(source.fetch_header());
        }
        
        auto calls = source.fetch_records();
        
        for (auto& call : calls) {
            auto reads = read_manager_.get().fetch_reads(mapped_region(call));
            
            const auto rmq = rmq_mapping_quality(reads);
            
            if (rmq < 20) {
                VcfRecord::Builder cb {call};
                
                cb.add_filter("MQ");
                
                call = cb.build();
            }
        }
        
        write(calls, dest);
    }
} // namespace Octopus
