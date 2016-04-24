//
//  vcf_record_factory.cpp
//  Octopus
//
//  Created by Daniel Cooke on 21/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "vcf_record_factory.hpp"

#include <string>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <iostream>

#include "genomic_region.hpp"
#include "mappable.hpp"
#include "mappable_algorithms.hpp"
#include "allele.hpp"
#include "variant_call.hpp"
#include "read_utils.hpp"
#include "string_utils.hpp"
#include "maths.hpp"

namespace Octopus
{
    VcfRecordFactory::VcfRecordFactory(const ReferenceGenome& reference, const ReadMap& reads,
                                       std::vector<SampleIdType> samples, bool sites_only)
    :
    reference_ {reference},
    reads_ {reads},
    samples_ {std::move(samples)},
    sites_only_ {sites_only}
    {}
    
//    auto count_alt_alleles(const Call& call)
//    {
//        std::unordered_map<Allele, unsigned> allele_counts {};
//        allele_counts.reserve(variants.size());
//        
//        for (const auto& allele : genotype_call.genotype) {
//            ++allele_counts[allele];
//        }
//        
//        std::vector<unsigned> result {};
//        result.reserve(variants.size());
//        
//        for (const auto& variant : variants) {
//            result.push_back(allele_counts[variant.get_alt_allele()]);
//        }
//        
//        return result;
//    }
//    
//    unsigned count_alleles(const Call& call)
//    {
//        std::unordered_set<Allele> unique_alleles {};
//        
//        for (const auto& allele : genotype_call.genotype) {
//            unique_alleles.emplace(allele);
//        }
//        
//        return static_cast<unsigned>(unique_alleles.size());
//    }
    
    void set_alt_alleles(const Call* call, VcfRecord::Builder& record)
    {
        if (auto p = dynamic_cast<const VariantCall*>(call)) {
            record.set_alt_allele(p->get_alternative().get_sequence());
        } else {
            record.set_refcall();
        }
    }
    
    void set_vcf_genotype(const SampleIdType& sample, const Call::GenotypeCall& genotype_call,
                          VcfRecord::Builder& record)
    {
        std::vector<VcfRecord::SequenceType> result {};
        result.reserve(genotype_call.genotype.ploidy());
        
        for (const auto& allele : genotype_call.genotype) {
            result.push_back(allele.get_sequence());
        }
        
        record.add_genotype(sample, result, VcfRecord::Builder::Phasing::Phased);
    }
    
    VcfRecord VcfRecordFactory::make(const std::unique_ptr<Call>& call) const
    {
        using std::to_string;
        
        auto result = VcfRecord::Builder {};
        
        const auto phred_quality = Maths::probability_to_phred<float>(call->get_quality(), 2);
        
        const auto& reference_allele = call->get_reference();
        
        const auto& region = call->get_region();
        
        result.set_chromosome(contig_name(region));
        result.set_position(region_begin(region));
        result.set_ref_allele(reference_allele.get_sequence());
        result.set_quality(phred_quality);
        
        set_alt_alleles(call.get(), result);
        
        //result.add_info("AC",  to_strings(count_alt_alleles(*call)));
        //result.add_info("AN",  to_string(count_alleles(*call)));
        
        result.add_info("NS",  to_string(count_samples_with_coverage(reads_, region)));
        result.add_info("DP",  to_string(sum_max_coverages(reads_, region)));
        result.add_info("SB",  Octopus::to_string(strand_bias(reads_, region), 2));
        result.add_info("BQ",  to_string(static_cast<unsigned>(rmq_base_quality(reads_, region))));
        result.add_info("MQ",  to_string(static_cast<unsigned>(rmq_mapping_quality(reads_, region))));
        result.add_info("MQ0", to_string(count_mapq_zero(reads_, region)));
        
        if (!sites_only_) {
            if (call->all_phased()) {
                result.set_format({"GT", "FT", "GQ", "DP", "BQ", "MQ", "PS", "PQ"});
            } else {
                result.set_format({"GT", "FT", "GQ", "DP", "BQ", "MQ"});
            }
            
            for (const auto& sample : samples_) {
                const auto& genotype_call = call->get_genotype_call(sample);
                
                set_vcf_genotype(sample, genotype_call, result);
                
                result.add_genotype_field(sample, "FT", "."); // TODO
                result.add_genotype_field(sample, "GQ", Octopus::to_string(Maths::probability_to_phred<float>(genotype_call.posterior), 2));
                result.add_genotype_field(sample, "DP", to_string(max_coverage(reads_.at(sample), region)));
                result.add_genotype_field(sample, "BQ", to_string(static_cast<unsigned>(rmq_base_quality(reads_.at(sample), region))));
                result.add_genotype_field(sample, "MQ", to_string(static_cast<unsigned>(rmq_mapping_quality(reads_.at(sample), region))));
                
                if (call->is_phased(sample)) {
                    const auto& phase = *genotype_call.phase;
                    result.add_genotype_field(sample, "PS", to_string(region_begin(phase.region) + 1));
                    result.add_genotype_field(sample, "PQ", Octopus::to_string(Maths::probability_to_phred<float>(phase.score), 2)); // TODO
                }
            }
        }
        
        return result.build_once();
    }
    
    namespace
    {
        struct CallWrapper : public Mappable<CallWrapper>
        {
            CallWrapper(const std::unique_ptr<Call>& call) : call {std::cref(call) } {}
            operator const std::unique_ptr<Call>&() const noexcept { return call.get(); }
            std::unique_ptr<Call>::pointer operator->() const noexcept { return call.get().get(); };
            std::reference_wrapper<const std::unique_ptr<Call>> call;
            const GenomicRegion& get_region() const noexcept { return call.get()->get_region(); }
        };
        
        auto wrap(const std::vector<std::unique_ptr<Call>>& calls)
        {
            return std::vector<CallWrapper> {std::begin(calls), std::end(calls)};
        }
    } // namespace
    
    template <typename ForwardIt>
    void parsimonise_mutually_exclusive_calls(ForwardIt first, ForwardIt last,
                                              std::vector<VcfRecord>& result)
    {
//        std::transform(first, last, std::ostream_iterator<std::string>(std::cout, "\n"),
//                       [] (const CallWrapper& call) {
//                           return to_string(call->get_region());
//                       });
    }
    
    std::vector<VcfRecord>
    VcfRecordFactory::make(const std::vector<std::unique_ptr<Call>>& calls) const
    {
        const auto wrapped_calls = wrap(calls);
        
        std::vector<VcfRecord> result {};
        result.reserve(calls.size());
        
        auto it = std::cbegin(wrapped_calls);
        
        while (it != std::cend(wrapped_calls)) {
            const auto it2 = find_first_overlapped(it, std::cend(wrapped_calls));
            
            std::transform(it, it2, std::back_inserter(result),
                           [this] (const auto& call) { return this->make(call); });
            
            if (it2 == std::cend(wrapped_calls)) break;
            
            it = find_first_not_overlapped(std::next(it2), std::cend(wrapped_calls), *it2);
            
            parsimonise_mutually_exclusive_calls(it2, it, result);
        }
        
        result.shrink_to_fit();
        
        return result;
    }
} // namespace Octopus
