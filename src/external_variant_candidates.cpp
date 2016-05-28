//
//  external_variant_candidates.cpp
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "external_variant_candidates.hpp"

#include <cstddef>
#include <algorithm>
#include <utility>

#include "vcf_record.hpp"
#include "variant.hpp"

namespace Octopus {

ExternalCandidateVariantGenerator::ExternalCandidateVariantGenerator(boost::filesystem::path path)
:
reader_ {std::make_shared<VcfReader>(std::move(path))}
{}

ExternalCandidateVariantGenerator::ExternalCandidateVariantGenerator(std::unique_ptr<VcfReader> reader)
:
reader_ {std::move(reader)}
{}

ExternalCandidateVariantGenerator::ExternalCandidateVariantGenerator(const std::shared_ptr<VcfReader>& reader)
:
reader_ {reader}
{}

std::vector<GenomicRegion> get_batch_regions(const GenomicRegion& region, VcfReader& reader,
                                             std::size_t max_batch_size)
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
    
    std::size_t max_batch_size {1000000000};
    
    // TODO: we really need iterators in VcfReader
    
    auto batches = get_batch_regions(region, reader, max_batch_size);
    
    for (const auto& batch : batches) {
        auto records = reader.fetch_records(batch, VcfReader::UnpackPolicy::Sites);
        
        for (const auto& record : records) {
            for (const auto& alt_allele : record.alt()) {
                const auto& ref_allele = record.ref();
                
                if (ref_allele.size() != alt_allele.size()) {
                    auto begin = record.pos();
                    
                    const auto p = std::mismatch(std::cbegin(ref_allele), std::cend(ref_allele),
                                                 std::cbegin(alt_allele), std::cend(alt_allele));
                    
                    Variant::SequenceType new_ref_allele {p.first, std::cend(ref_allele)};
                    Variant::SequenceType new_alt_allele {p.second, std::cend(alt_allele)};
                    
                    begin += std::distance(std::cbegin(ref_allele), p.first);
                    
                    result.emplace_back(record.chrom(), begin - 1,
                                        std::move(new_ref_allele), std::move(new_alt_allele));
                } else {
                    result.emplace_back(record.chrom(), record.pos() - 1, record.ref(), alt_allele);
                }
            }
        }
    }
    
    result.shrink_to_fit();
    
    std::sort(std::begin(result), std::end(result));
    
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    
    return result;
}

std::vector<Variant> ExternalCandidateVariantGenerator::generate_candidates(const GenomicRegion& region)
{
    return fetch_variants(region, *reader_);
}

} // namespace Octopus
