//
//  variant_call_filter.cpp
//  Octopus
//
//  Created by Daniel Cooke on 31/05/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "variant_call_filter.hpp"

#include <unordered_map>

#include "common.hpp"

#include "vcf_reader.hpp"
#include "vcf_writer.hpp"
#include "vcf_header.hpp"
#include "vcf_record.hpp"

#include "genomic_region.hpp"
#include "mappable_flat_set.hpp"
#include "mappable_map.hpp"
#include "aligned_read.hpp"
#include "read_utils.hpp"

#include "genotype_reader.hpp"
#include "haplotype_likelihood_cache.hpp"

namespace Octopus
{
VariantCallFilter::VariantCallFilter(const ReferenceGenome& reference, const ReadManager& read_manager)
:
reference_ {reference},
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

using PhaseMap = MappableMap<ContigNameType, GenomicRegion, MappableFlatSet<GenomicRegion>>;

using SamplePhaseMap = std::unordered_map<SampleIdType, PhaseMap>;

SamplePhaseMap extract_phase_map(const VcfReader& calls)
{
    const auto samples = calls.fetch_header().samples();
    
    SamplePhaseMap result {};
    result.reserve(samples.size());
    
    
    
    return result;
}

void VariantCallFilter::filter(const VcfReader& source, VcfWriter& dest, const InputRegionMap& regions)
{
    if (!dest.is_header_written()) {
        dest << source.fetch_header();
    }
    
    GenotypeReader gr {reference_, source};
    
    auto genotypes = gr.extract_genotype(reference_.get().contig_region("22"));
    
    for (const auto& p : genotypes) {
        for (const auto& g : p.second) {
            std::cout << g << std::endl;
        }
    }
    
    for (const auto& p : regions) {
        for (const auto& region : p.second) {
            auto calls = source.fetch_records(region);
            
            for (auto& call : calls) {
                const auto call_region = mapped_region(call);
                
                auto reads = read_manager_.get().fetch_reads(call_region);
                
                VcfRecord::Builder cb {call};
                
                bool filtered {false};
                
                if (call.qual() && *call.qual() < 10) {
                    cb.add_filter("q10");
                    filtered = true;
                }
                
                if (rmq_mapping_quality(reads) < 40) {
                    cb.add_filter("MQ");
                    filtered = true;
                }
                
                if (!filtered) {
                    cb.set_passed();
                }
                
                call = cb.build();
            }
            
            write(calls, dest);
        }
    }
}
} // namespace Octopus
