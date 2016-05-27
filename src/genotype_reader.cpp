//
//  genotype_reader.cpp
//  Octopus
//
//  Created by Daniel Cooke on 23/05/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "genotype_reader.hpp"

#include <utility>
#include <functional>
#include <iterator>
#include <algorithm>
#include <cassert>

#include <boost/lexical_cast.hpp>

#include "vcf_header.hpp"
#include "vcf_record.hpp"
#include "contig_region.hpp"
#include "mappable.hpp"
#include "mappable_algorithms.hpp"
#include "allele.hpp"

#include <iostream> // DEBUG

namespace Octopus
{
GenotypeReader::GenotypeReader(const ReferenceGenome& reference, VcfReader&& variant_reader)
:
reference_ {reference},
variant_reader_ {std::move(variant_reader)},
samples_ {variant_reader_.fetch_header().samples()}
{}

namespace
{
    auto extract_contig_region(const VcfRecord& call)
    {
        const auto begin = static_cast<ContigRegion::SizeType>(call.position()) - 1;
        return ContigRegion {
            begin, begin + static_cast<ContigRegion::SizeType>(call.ref_allele().size())
        };
    }
    
    auto extract_phase_region(const VcfRecord& call, const SampleIdType& sample)
    {
        if (call.is_sample_phased(sample)) {
            assert(call.has_format("PS"));
            return GenomicRegion {
                call.chromosome_name(),
                boost::lexical_cast<ContigRegion::SizeType>(call.get_sample_value(sample, "PS").front()) - 1,
                static_cast<ContigRegion::SizeType>(call.position() + call.ref_allele().size()) - 1
            };
        }
        return GenomicRegion {call.chromosome_name(), extract_contig_region(call)};
    }
    
    struct CallWrapper : public Mappable<CallWrapper>
    {
        CallWrapper(const VcfRecord& record, const SampleIdType& sample)
        :
        call {std::cref(record)}, phase_region {extract_phase_region(record, sample)} {}
        
        std::reference_wrapper<const VcfRecord> call;
        GenomicRegion phase_region;
        const GenomicRegion& mapped_region() const noexcept { return phase_region; }
    };
    
    auto wrap_calls(const std::vector<VcfRecord>& calls, const SampleIdType& sample)
    {
        std::vector<CallWrapper> result {};
        result.reserve(calls.size());
        
        for (const auto& call : calls) {
            result.emplace_back(call, sample);
        }
        
        return result;
    }
    
    decltype(auto) extract_genotype(const CallWrapper& call, const SampleIdType& sample)
    {
        return call.call.get().get_sample_value(sample, "GT");
    }
    
    auto extract_ploidy(const std::vector<CallWrapper>& phased_calls,
                        const SampleIdType& sample)
    {
        assert(!phased_calls.empty());
        return extract_genotype(phased_calls.front(), sample).size();
    }
    
    auto make_genotype(std::vector<Haplotype::Builder>&& haplotypes)
    {
        Genotype<Haplotype> result {static_cast<unsigned>(haplotypes.size())};
        
        for (auto& haplotype : haplotypes) {
            result.emplace(haplotype.build());
        }
        
        return result;
    }
    
    auto extract_contig_region(const CallWrapper& call)
    {
        return extract_contig_region(call.call.get());
    }
    
    Genotype<Haplotype> extract_genotype(const std::vector<CallWrapper>& phased_calls,
                                         const SampleIdType& sample, const ReferenceGenome& reference)
    {
        assert(!phased_calls.empty());
        
        const auto ploidy = extract_ploidy(phased_calls, sample);
        
        const auto region = encompassing_region(phased_calls);
        
        std::vector<Haplotype::Builder> haplotypes(ploidy, Haplotype::Builder {region, reference});
        
        for (const auto& call : phased_calls) {
            const auto& genotype = extract_genotype(call, sample);
            
            for (unsigned i {0}; i < ploidy; ++i) {
                haplotypes[i].push_back(ContigAllele {extract_contig_region(call), genotype[i]});
            }
        }
        
        return make_genotype(std::move(haplotypes));
    }
} // namespace

GenotypeReader::GenotypeMap GenotypeReader::extract_genotype(const GenomicRegion& region)
{
    const auto calls = variant_reader_.fetch_records(region);
    
    GenotypeMap result {samples_.size()};
    
    for (const auto& sample : samples_) {
        const auto wrapped_calls = wrap_calls(calls, sample);
        
        for (const auto& phased_calls : segment_overlapped_copy(wrapped_calls)) {
            result[sample].insert(::Octopus::extract_genotype(phased_calls, sample, reference_.get()));
        }
    }
    
    return result;
}
} // namespace Octopus
