//
//  external_variant_candidates.cpp
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "external_variant_candidates.hpp"

#include <cstddef>

#include "vcf_record.hpp"
#include "variant.hpp"

ExternalVariantCandidates::ExternalVariantCandidates(VcfReader& reader)
:
reader_ {reader}
{}

std::vector<GenomicRegion> get_batch_regions(const GenomicRegion& region, const VcfReader& reader, std::size_t max_batch_size)
{
    std::vector<GenomicRegion> result {};
    
    if (reader.count_records(region) > max_batch_size) {
        // umm?
    } else {
        result.push_back(region);
    }
    
    return result;
}

std::vector<Variant> fetch_variants(const GenomicRegion& region, VcfReader& reader)
{
    std::vector<Variant> result {};
    result.reserve(reader.count_records(region));
    
    size_t max_batch_size {10000};
    
    auto batches = get_batch_regions(region, reader, max_batch_size);
    
    for (const auto& batch : batches) {
        auto vcf_records = reader.fetch_records(batch, VcfReader::Unpack::AllButSamples);
        for (const auto& record : vcf_records) {
            for (const auto& alt_allele : record.get_alt_alleles()) {
                result.emplace_back(record.get_chromosome_name(), record.get_position(), record.get_ref_allele(), alt_allele);
            }
        }
    }
    
    return result;
}

std::vector<Variant> ExternalVariantCandidates::get_candidates(const GenomicRegion& region)
{
    return fetch_variants(region, reader_);
}
