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

#include <boost/optional.hpp>

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
    
    void set_alt_allele(const Call* call, VcfRecord::Builder& record,
                        boost::optional<char> first_reference_base)
    {
        if (auto p = dynamic_cast<const VariantCall*>(call)) {
            auto sequence = p->get_alternative().get_sequence();
            if (first_reference_base) sequence.front() = *first_reference_base;
            record.set_alt_allele(std::move(sequence));
        } else {
            record.set_refcall();
        }
    }
    
    void set_vcf_genotype(const SampleIdType& sample, const Call::GenotypeCall& genotype_call,
                          VcfRecord::Builder& record, boost::optional<char> first_reference_base)
    {
        std::vector<VcfRecord::SequenceType> result {};
        result.reserve(genotype_call.genotype.ploidy());
        
        for (const auto& allele : genotype_call.genotype) {
            auto sequence = allele.get_sequence();
            if (first_reference_base) sequence.front() = *first_reference_base;
            result.push_back(std::move(sequence));
        }
        
        record.add_genotype(sample, result, VcfRecord::Builder::Phasing::Phased);
    }
    
    VcfRecord VcfRecordFactory::make(std::unique_ptr<Call> call) const
    {
        using std::to_string;
        
        auto result = VcfRecord::Builder {};
        
        const auto phred_quality = Maths::probability_to_phred<float>(call->get_quality(), 2);
        
        const auto& region = call->get_region();
        
        auto reference_allele = call->get_reference().get_sequence();
        
        boost::optional<char> first_reference_base {};
        
        if (reference_allele.front() == '#') {
            first_reference_base = reference_.get_sequence(head_position(region)).front();
            reference_allele.front() = *first_reference_base;
        }
        
        result.set_chromosome(contig_name(region));
        result.set_position(region_begin(region));
        result.set_ref_allele(std::move(reference_allele));
        result.set_quality(phred_quality);
        
        set_alt_allele(call.get(), result, first_reference_base);
        
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
                
                set_vcf_genotype(sample, genotype_call, result, first_reference_base);
                
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
            CallWrapper(std::unique_ptr<Call> call) : call {std::move(call) } {}
            
            CallWrapper(const CallWrapper&)            = delete;
            CallWrapper& operator=(const CallWrapper&) = delete;
            CallWrapper(CallWrapper&&)                 = default;
            CallWrapper& operator=(CallWrapper&&)      = default;
            
            operator const std::unique_ptr<Call>&() const noexcept { return call; }
            operator std::unique_ptr<Call>&() noexcept { return call; }
            
            std::unique_ptr<Call>::pointer operator->() const noexcept { return call.get(); };
            
            const GenomicRegion& get_region() const noexcept { return call->get_region(); }
            
            std::unique_ptr<Call> call;
        };
        
        template <typename T>
        auto wrap(std::vector<std::unique_ptr<T>>&& calls)
        {
            return std::vector<CallWrapper> {
                std::make_move_iterator(std::begin(calls)),
                std::make_move_iterator(std::end(calls))
            };
        }
    } // namespace
    
    bool are_in_phase(const Call::GenotypeCall& lhs, const Call::GenotypeCall& rhs)
    {
        return overlaps(lhs.phase->region, rhs.genotype);
    }
    
    std::vector<VcfRecord>
    VcfRecordFactory::make(std::vector<std::unique_ptr<Call>>&& calls) const
    {
        auto wrapped_calls = wrap(std::move(calls));
        
        for (auto& call : wrapped_calls) {
            call->parsimonise('#');
        }
        
        std::vector<VcfRecord> result {};
        result.reserve(calls.size());
        
        auto it = std::begin(wrapped_calls);
        
        while (it != std::end(wrapped_calls)) {
            const auto it2 = find_first_overlapped(it, std::end(wrapped_calls));
            
            std::transform(std::make_move_iterator(it), std::make_move_iterator(it2),
                           std::back_inserter(result),
                           [this] (CallWrapper&& call) {
                               return this->make(std::move(call.call));
                           });
            
            if (it2 == std::end(wrapped_calls)) break;
            
            it = find_first_not_overlapped(std::next(it2), std::end(wrapped_calls), *it2);
            
            if ((*it2)->get_reference().get_sequence().front() == '#') {
                for (const auto& sample : samples_) {
                    const auto& genotype_call = (*it2)->get_genotype_call(sample);
                    
                    const auto& old_genotype = genotype_call.genotype;
                    
                    const auto ploidy = old_genotype.ploidy();
                    
                    const auto actual_reference_base = reference_.get_sequence(head_position(*it2)).front();
                    
                    Genotype<Allele> new_genotype {ploidy};
                    
                    for (unsigned i {0}; i < ploidy; ++i) {
                        if (old_genotype[i].get_sequence().front() == '#') {
                            auto new_sequence = old_genotype[i].get_sequence();
                            new_sequence.front() = actual_reference_base;
                            Allele new_allele {(*it2).get_region(), std::move(new_sequence)};
                            new_genotype.emplace(std::move(new_allele));
                        } else {
                            new_genotype.emplace(old_genotype[i]);
                        }
                    }
                    
                    std::cout << old_genotype << " has been replaced with " << new_genotype << std::endl;
                }
            }
            
            auto rit1 = std::make_reverse_iterator(it);
            const auto rit2 = std::prev(std::make_reverse_iterator(it2));
            
            for (; rit1 != rit2; ++rit1) {
                const auto& curr_call = *rit1;
                
                bool has_missing_allele {false};
                
                if (curr_call->get_reference().get_sequence().front() == '#') {
                    const auto actual_reference_base = reference_.get_sequence(head_position(curr_call)).front();
                    
                    const auto& prev_call = *std::next(rit1);
                    
                    if (overlaps(curr_call, prev_call)) {
                        for (const auto& sample : samples_) {
                            const auto& genotype_call = curr_call->get_genotype_call(sample);
                            
                            const auto& old_genotype = genotype_call.genotype;
                            
                            const auto& prev_genotype_call = prev_call->get_genotype_call(sample);
                            const auto& prev_genotype = prev_genotype_call.genotype;
                            
                            const auto ploidy = old_genotype.ploidy();
                            
                            Genotype<Allele> new_genotype {ploidy};
                            
                            if (are_in_phase(prev_genotype_call, genotype_call)) {
                                for (unsigned i {0}; i < ploidy; ++i) {
                                    if (old_genotype[i].get_sequence().front() == '#') {
                                        if (splice(prev_genotype[i], head_position(curr_call)).get_sequence().front() == actual_reference_base) {
                                            new_genotype.emplace(old_genotype[i]);
                                        } else {
                                            auto new_sequence = old_genotype[i].get_sequence();
                                            new_sequence.front() = '*';
                                            Allele new_allele {curr_call.get_region(), std::move(new_sequence)};
                                            new_genotype.emplace(std::move(new_allele));
                                            has_missing_allele = true;
                                        }
                                    } else {
                                        new_genotype.emplace(old_genotype[i]);
                                    }
                                }
                            } else {
                                has_missing_allele = true;
                                for (unsigned i {0}; i < ploidy; ++i) {
                                    if (old_genotype[i].get_sequence().front() == '#') {
                                        auto new_sequence = old_genotype[i].get_sequence();
                                        new_sequence.front() = '*';
                                        Allele new_allele {curr_call.get_region(), std::move(new_sequence)};
                                        new_genotype.emplace(std::move(new_allele));
                                    } else {
                                        new_genotype.emplace(old_genotype[i]);
                                    }
                                }
                            }
                            
                            std::cout << old_genotype << " has been replaced with " << new_genotype << std::endl;
                        }
                    } else {
                        for (const auto& sample : samples_) {
                            const auto& genotype_call = curr_call->get_genotype_call(sample);
                            
                            const auto& old_genotype = genotype_call.genotype;
                            
                            const auto ploidy = old_genotype.ploidy();
                            
                            Genotype<Allele> new_genotype {ploidy};
                            
                             for (unsigned i {0}; i < ploidy; ++i) {
                                 if (old_genotype[i].get_sequence().front() == '#') {
                                     auto new_sequence = old_genotype[i].get_sequence();
                                     new_sequence.front() = actual_reference_base;
                                     Allele new_allele {curr_call.get_region(), std::move(new_sequence)};
                                     new_genotype.emplace(std::move(new_allele));
                                 } else {
                                     new_genotype.emplace(old_genotype[i]);
                                 }
                             }
                            
                            std::cout << old_genotype << " has been replaced with " << new_genotype << std::endl;
                        }
                    }
                    
                    
                }
            }
            
            auto segements = segment_by_begin_copy(std::make_move_iterator(it2),
                                                   std::make_move_iterator(it));
            
            
            
//            std::transform(std::make_move_iterator(it2), std::make_move_iterator(it),
//                           std::back_inserter(result),
//                           [this] (CallWrapper&& call) {
//                               return this->make(std::move(call.call));
//                           });
        }
        
        result.shrink_to_fit();
        
        return result;
    }
    
    // private methods
    
//    VcfRecord VcfRecordFactory::make_single(const std::vector<std::unique_ptr<Call>>& calls) const
//    {
//        return VcfRecord {};
//    }
} // namespace Octopus
