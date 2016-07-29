//
//  source_file_reader.cpp
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "source_file_reader.hpp"

#include <cstddef>
#include <algorithm>
#include <utility>

#include "vcf_record.hpp"

namespace octopus { namespace core { namespace generators
{
SourceFileReader::SourceFileReader(boost::filesystem::path path)
:
reader_ {std::make_shared<VcfReader>(std::move(path))}
{}

SourceFileReader::SourceFileReader(std::unique_ptr<const VcfReader> reader)
:
reader_ {std::move(reader)}
{}

SourceFileReader::SourceFileReader(const std::shared_ptr<const VcfReader>& reader)
:
reader_ {reader}
{}

std::vector<GenomicRegion> get_batch_regions(const GenomicRegion& region, const VcfReader& reader,
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

static bool is_missing(const VcfRecord::NucleotideSequence& allele)
{
    return allele == "*";
}

std::vector<Variant> fetch_variants(const GenomicRegion& region, const VcfReader& reader)
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
                if (!is_missing(alt_allele)) {
                    const auto& ref_allele = record.ref();
                    
                    if (ref_allele.size() != alt_allele.size()) {
                        auto begin = record.pos();
                        
                        const auto p = std::mismatch(std::cbegin(ref_allele), std::cend(ref_allele),
                                                     std::cbegin(alt_allele), std::cend(alt_allele));
                        
                        Variant::NucleotideSequence new_ref_allele {p.first, std::cend(ref_allele)};
                        Variant::NucleotideSequence new_alt_allele {p.second, std::cend(alt_allele)};
                        
                        begin += std::distance(std::cbegin(ref_allele), p.first);
                        
                        result.emplace_back(record.chrom(), begin - 1,
                                            std::move(new_ref_allele), std::move(new_alt_allele));
                    } else {
                        result.emplace_back(record.chrom(), record.pos() - 1, record.ref(), alt_allele);
                    }
                }
            }
        }
    }
    
    result.shrink_to_fit();
    
    std::sort(std::begin(result), std::end(result));
    
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    
    return result;
}

std::vector<Variant> SourceFileReader::generate_candidates(const GenomicRegion& region)
{
    return fetch_variants(region, *reader_);
}
} // namespace generators
} // namespace core
} // namespace octopus
