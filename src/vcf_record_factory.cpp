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
        
        for (auto it = std::begin(wrapped_calls); it != std::end(wrapped_calls);) {
            if (is_empty(it->get_region())) {
                auto it2 = std::find_if_not(std::next(it), std::end(wrapped_calls),
                                            [it] (const auto& call) {
                                                return call->get_region() == it->get_region();
                                            });
                
                if (it2 == std::end(wrapped_calls)) break;
                
                if (!overlaps(*it, *it2)) {
                    it = it2;
                    continue;
                }
                
                auto it3 = find_first_not_overlapped(std::next(it2), std::end(wrapped_calls), *it);
                
                // now everything between it and it2 is an insertion, anything between
                // it2 and it3 is another call which will have inserted sequence we want to remove.
                // Note the genotype calls of all insertions must be the same as they are in the
                // same region
                
                std::for_each(it2, it3, [this, it] (auto& call) {
                    for (const auto& sample : samples_) {
                        const auto& insertion_genotype = (*it)->get_genotype_call(sample).genotype;
                        
                        auto& sample_genotype = call->get_genotype_call(sample).genotype;
                        
                        std::vector<Allele::SequenceType> resolved_alleles {};
                        resolved_alleles.reserve(sample_genotype.ploidy());
                        
                        std::transform(std::cbegin(sample_genotype), std::cend(sample_genotype),
                                       std::cbegin(insertion_genotype), std::back_inserter(resolved_alleles),
                                       [] (const Allele& allele1, const Allele& allele2) {
                                           if (is_insertion(allele2)) {
                                               const auto& old_sequence = allele1.get_sequence();
                                               return Allele::SequenceType {
                                                   std::next(std::cbegin(old_sequence), sequence_size(allele2)),
                                                   std::cend(old_sequence)
                                               };
                                           }
                                           return allele1.get_sequence();
                                       });
                        
                        Genotype<Allele> new_genotype {sample_genotype.ploidy()};
                        
                        for (auto& sequence : resolved_alleles) {
                            new_genotype.emplace(Allele {mapped_region(sample_genotype), std::move(sequence)});
                        }
                        
                        sample_genotype = std::move(new_genotype);
                    }
                });
                
                it = it3;
            } else {
                ++it;
            }
        }
        
        for (auto& call : wrapped_calls) {
            call->parsimonise('#');
        }
        
        // TODO: we need to adjust the phase regions of all calls where the modified parsimonise
        // call was the the leftmost call within the phase region, otherwise the phase sets
        // won't line up properly.
        
        std::vector<VcfRecord> result {};
        result.reserve(calls.size());
        
        for (auto it = std::begin(wrapped_calls); it != std::end(wrapped_calls);) {
            const auto it2 = find_first_overlapped(it, std::end(wrapped_calls));
            
            std::transform(std::make_move_iterator(it), std::make_move_iterator(it2),
                           std::back_inserter(result),
                           [this] (CallWrapper&& call) {
                               call->replace('#', reference_.get_sequence(head_position(call->get_region())).front());
                               return this->make(std::move(call.call));
                           });
            
            if (it2 == std::end(wrapped_calls)) break;
            
            it = find_first_not_overlapped(std::next(it2), std::end(wrapped_calls), *it2);
            
            if ((*it2)->get_reference().get_sequence().front() == '#') {
                const auto actual_reference_base = reference_.get_sequence(head_position(*it2)).front();
                
                auto new_sequence = (*it2)->get_reference().get_sequence();
                new_sequence.front() = actual_reference_base;
                Allele new_allele {(*it2).get_region(), std::move(new_sequence)};
                
                (*it2)->replace((*it2)->get_reference(), std::move(new_allele));
                
                std::unordered_map<Allele, Allele> replacements {};
                
                for (const auto& sample : samples_) {
                    auto& genotype_call = (*it2)->get_genotype_call(sample);
                    
                    auto& old_genotype = genotype_call.genotype;
                    
                    const auto ploidy = old_genotype.ploidy();
                    
                    Genotype<Allele> new_genotype {ploidy};
                    
                    for (unsigned i {0}; i < ploidy; ++i) {
                        if (old_genotype[i].get_sequence().front() == '#') {
                            auto new_sequence = old_genotype[i].get_sequence();
                            new_sequence.front() = actual_reference_base;
                            Allele new_allele {(*it2).get_region(), std::move(new_sequence)};
                            replacements.emplace(old_genotype[i], new_allele);
                            new_genotype.emplace(std::move(new_allele));
                        } else {
                            new_genotype.emplace(old_genotype[i]);
                        }
                    }
                    
                    old_genotype = std::move(new_genotype);
                }
                
                for (auto& p : replacements) {
                    (*it2)->replace(p.first, p.second);
                }
            }
            
            for (auto it3 = std::next(it2); it3 != it; ++it3) {
                auto& curr_call = *it3;
                
                bool has_missing_allele {false};
                
                std::unordered_map<Allele, Allele> replacements {};
                
                if (curr_call->get_reference().get_sequence().front() == '#') {
                    const auto actual_reference_base = reference_.get_sequence(head_position(curr_call)).front();
                    
                    auto new_sequence = curr_call->get_reference().get_sequence();
                    new_sequence.front() = actual_reference_base;
                    Allele new_allele {curr_call.get_region(), std::move(new_sequence)};
                    
                    curr_call->replace(curr_call->get_reference(), std::move(new_allele));
                    
                    const auto& prev_call = *std::prev(it3);
                    
                    if (overlaps(curr_call, prev_call)) {
                        for (const auto& sample : samples_) {
                            auto& genotype_call = curr_call->get_genotype_call(sample);
                            
                            auto& old_genotype = genotype_call.genotype;
                            
                            const auto& prev_genotype_call = prev_call->get_genotype_call(sample);
                            const auto& prev_genotype = prev_genotype_call.genotype;
                            
                            const auto ploidy = old_genotype.ploidy();
                            
                            Genotype<Allele> new_genotype {ploidy};
                            
                            if (are_in_phase(prev_genotype_call, genotype_call)) {
                                for (unsigned i {0}; i < ploidy; ++i) {
                                    if (old_genotype[i].get_sequence().empty()) {
                                        Allele::SequenceType new_sequence(size(curr_call.get_region()), '*');
                                        Allele new_allele {curr_call.get_region(), std::move(new_sequence)};
                                        new_genotype.emplace(std::move(new_allele));
                                    } else if (old_genotype[i].get_sequence().front() == '#') {
                                        if (splice(prev_genotype[i], head_position(curr_call)).get_sequence().front() == actual_reference_base) {
                                            auto new_sequence = old_genotype[i].get_sequence();
                                            new_sequence.front() = actual_reference_base;
                                            Allele new_allele {(*it2).get_region(), std::move(new_sequence)};
                                            replacements.emplace(old_genotype[i], new_allele);
                                            new_genotype.emplace(std::move(new_allele));
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
                                    if (old_genotype[i].get_sequence().empty()) {
                                        Allele::SequenceType new_sequence(size(curr_call.get_region()), '*');
                                        Allele new_allele {curr_call.get_region(), std::move(new_sequence)};
                                        new_genotype.emplace(std::move(new_allele));
                                    } else if (old_genotype[i].get_sequence().front() == '#') {
                                        auto new_sequence = old_genotype[i].get_sequence();
                                        new_sequence.front() = '.';
                                        Allele new_allele {curr_call.get_region(), std::move(new_sequence)};
                                        new_genotype.emplace(std::move(new_allele));
                                    } else {
                                        new_genotype.emplace(old_genotype[i]);
                                    }
                                }
                            }
                            
                            old_genotype = std::move(new_genotype);
                        }
                    } else {
                        for (const auto& sample : samples_) {
                            auto& genotype_call = curr_call->get_genotype_call(sample);
                            
                            auto& old_genotype = genotype_call.genotype;
                            
                            const auto ploidy = old_genotype.ploidy();
                            
                            Genotype<Allele> new_genotype {ploidy};
                            
                            for (unsigned i {0}; i < ploidy; ++i) {
                                if (old_genotype[i].get_sequence().front() == '#') {
                                    auto new_sequence = old_genotype[i].get_sequence();
                                    new_sequence.front() = actual_reference_base;
                                    Allele new_allele {curr_call.get_region(), std::move(new_sequence)};
                                    replacements.emplace(old_genotype[i], new_allele);
                                    new_genotype.emplace(std::move(new_allele));
                                } else {
                                    new_genotype.emplace(old_genotype[i]);
                                }
                            }
                            
                            old_genotype = std::move(new_genotype);
                        }
                    }
                } else {
                    for (const auto& sample : samples_) {
                        auto& genotype_call = curr_call->get_genotype_call(sample);
                        
                        auto& old_genotype = genotype_call.genotype;
                        
                        const auto ploidy = old_genotype.ploidy();
                        
                        Genotype<Allele> new_genotype {ploidy};
                        
                        for (unsigned i {0}; i < ploidy; ++i) {
                            if (old_genotype[i].get_sequence().empty()) {
                                Allele::SequenceType new_sequence(size(curr_call.get_region()), '*');
                                Allele new_allele {curr_call.get_region(), std::move(new_sequence)};
                                new_genotype.emplace(std::move(new_allele));
                            } else {
                                new_genotype.emplace(old_genotype[i]);
                            }
                        }
                        
                        old_genotype = std::move(new_genotype);
                    }
                }
                
                for (auto& p : replacements) {
                    curr_call->replace(p.first, p.second);
                }
            }
            
            std::for_each(it2, it, [] (auto& call) {
                call->replace_uncalled_genotype_alleles(Allele {call->get_region(), "."}, '*');
            });
            
            // At this point, all genotypes field contain normal bases, or '.' or '*', but not
            // '#'.
            
            auto segements = segment_by_begin_copy(std::make_move_iterator(it2),
                                                   std::make_move_iterator(it));
            
            for (auto&& segment : segements) {
                for (auto&& new_segment : segment_by_end_move(segment)) {
                    std::vector<std::unique_ptr<Call>> final_segment {};
                    std::transform(std::make_move_iterator(std::begin(new_segment)),
                                   std::make_move_iterator(std::end(new_segment)),
                                   std::back_inserter(final_segment),
                                   [] (auto&& call) -> std::unique_ptr<Call>&& {
                                       return std::move(call.call);
                                   });
                    result.emplace_back(this->make_segment(std::move(final_segment)));
                }
            }
        }
        
        result.shrink_to_fit();
        
        return result;
    }
    
    // private methods
    
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
    
    std::vector<VcfRecord::SequenceType>
    extract_all_genotyped_alleles(const Call* call, const std::vector<SampleIdType>& samples)
    {
        std::vector<VcfRecord::SequenceType> result {};
        
        for (const auto& sample : samples) {
            const auto& called_genotype = call->get_genotype_call(sample).genotype;
            
            std::transform(std::cbegin(called_genotype), std::cend(called_genotype),
                           std::back_inserter(result), [] (const Allele& allele) {
                               auto result = allele.get_sequence();
                               std::replace(std::begin(result), std::end(result), '*', '~');
                               return result;
                           });
        }
        
        const VcfRecord::SequenceType missing_allele(1, '.');
        
        auto it = std::remove(std::begin(result), std::end(result), missing_allele);
        
        result.erase(it, std::end(result));
        
        std::sort(std::begin(result), std::end(result));
        
        it = std::unique(std::begin(result), std::end(result));
        
        result.erase(it, std::end(result));
        
        for (auto& alt : result) {
            std::replace(std::begin(alt), std::end(alt), '~', '*');
        }
        
        return result;
    }
    
    void set_alt_alleles(const Call* call, VcfRecord::Builder& record,
                         const std::vector<SampleIdType>& samples)
    {
        auto alts = extract_all_genotyped_alleles(call, samples);
        
        auto it = std::find(std::begin(alts), std::end(alts), call->get_reference().get_sequence());
        
        if (it != std::end(alts)) {
            alts.erase(it);
        }
        
        record.set_alt_alleles(std::move(alts));
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
    
    VcfRecord VcfRecordFactory::make(std::unique_ptr<Call> call) const
    {
        using std::to_string;
        
        auto result = VcfRecord::Builder {};
        
        const auto phred_quality = Maths::probability_to_phred<float>(call->get_quality(), 2);
        
        const auto& region = call->get_region();
        
        result.set_chromosome(contig_name(region));
        result.set_position(region_begin(region));
        result.set_ref_allele(call->get_reference().get_sequence());
        result.set_quality(phred_quality);
        
        set_alt_alleles(call.get(), result, samples_);
        
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
    
    VcfRecord VcfRecordFactory::make_segment(std::vector<std::unique_ptr<Call>>&& calls) const
    {
        assert(!calls.empty());
        
        if (calls.size() == 1) {
            return make(std::move(calls.front()));
        }
        
        using std::to_string;
        
        auto result = VcfRecord::Builder {};
        
        const auto phred_quality = Maths::probability_to_phred<float>(calls.front()->get_quality(), 2);
        
        const auto& region = calls.front()->get_region();
        
        result.set_chromosome(contig_name(region));
        result.set_position(region_begin(region));
        result.set_ref_allele(calls.front()->get_reference().get_sequence());
        result.set_quality(phred_quality);
        
        std::vector<std::vector<VcfRecord::SequenceType>> resolved_genotypes {};
        resolved_genotypes.reserve(samples_.size());
        
        for (const auto& sample : samples_) {
            const auto ploidy = calls.front()->get_genotype_call(sample).genotype.ploidy();
            
            std::vector<VcfRecord::SequenceType> resolved_sample_genotype(ploidy);
            
            const auto& first_called_genotype = calls.front()->get_genotype_call(sample).genotype;
            
            std::transform(std::cbegin(first_called_genotype), std::cend(first_called_genotype),
                           std::begin(resolved_sample_genotype),
                           [] (const Allele& allele) {
                               return allele.get_sequence();
                           });
            
            std::for_each(std::next(std::cbegin(calls)), std::cend(calls), [&] (const auto& call) {
                const auto& called_genotype = call->get_genotype_call(sample).genotype;
                
                std::transform(std::cbegin(called_genotype), std::cend(called_genotype),
                               std::cbegin(resolved_sample_genotype),
                               std::begin(resolved_sample_genotype),
                               [] (const Allele& allele, const auto& curr) {
                                   const auto& seq = allele.get_sequence();
                                   
                                   if (seq.front() == '.' || seq.front() == '*') {
                                       return curr;
                                   }
                                   
                                   return seq;
                               });
            });
            
            resolved_genotypes.push_back(std::move(resolved_sample_genotype));
        }
        
        std::vector<VcfRecord::SequenceType> alt_alleles {};
        
        for (const auto& genotype : resolved_genotypes) {
            alt_alleles.insert(std::end(alt_alleles), std::cbegin(genotype), std::cend(genotype));
        }
        
        const VcfRecord::SequenceType missing_allele(1, '.');
        
        auto it = std::remove(std::begin(alt_alleles), std::end(alt_alleles), missing_allele);
        
        it = std::remove(std::begin(alt_alleles), it, calls.front()->get_reference().get_sequence());
        
        alt_alleles.erase(it, std::end(alt_alleles));
        
        std::sort(std::begin(alt_alleles), std::end(alt_alleles));
        
        it = std::unique(std::begin(alt_alleles), std::end(alt_alleles));
        
        alt_alleles.erase(it, std::end(alt_alleles));
        
        result.set_alt_alleles(std::move(alt_alleles));
        
        //result.add_info("AC",  to_strings(count_alt_alleles(*call)));
        //result.add_info("AN",  to_string(count_alleles(*call)));
        
        result.add_info("NS",  to_string(count_samples_with_coverage(reads_, region)));
        result.add_info("DP",  to_string(sum_max_coverages(reads_, region)));
        result.add_info("SB",  Octopus::to_string(strand_bias(reads_, region), 2));
        result.add_info("BQ",  to_string(static_cast<unsigned>(rmq_base_quality(reads_, region))));
        result.add_info("MQ",  to_string(static_cast<unsigned>(rmq_mapping_quality(reads_, region))));
        result.add_info("MQ0", to_string(count_mapq_zero(reads_, region)));
        
        if (!sites_only_) {
            if (calls.front()->all_phased()) {
                result.set_format({"GT", "FT", "GQ", "DP", "BQ", "MQ", "PS", "PQ"});
            } else {
                result.set_format({"GT", "FT", "GQ", "DP", "BQ", "MQ"});
            }
            
            auto sample_itr = std::begin(resolved_genotypes);
            
            for (const auto& sample : samples_) {
                const auto posterior = calls.front()->get_genotype_call(sample).posterior;
                
                result.add_genotype(sample, *sample_itr++, VcfRecord::Builder::Phasing::Phased);
                
                result.add_genotype_field(sample, "FT", "."); // TODO
                result.add_genotype_field(sample, "GQ", Octopus::to_string(Maths::probability_to_phred<float>(posterior), 2));
                result.add_genotype_field(sample, "DP", to_string(max_coverage(reads_.at(sample), region)));
                result.add_genotype_field(sample, "BQ", to_string(static_cast<unsigned>(rmq_base_quality(reads_.at(sample), region))));
                result.add_genotype_field(sample, "MQ", to_string(static_cast<unsigned>(rmq_mapping_quality(reads_.at(sample), region))));
                
                if (calls.front()->is_phased(sample)) {
                    const auto phase = *calls.front()->get_genotype_call(sample).phase;
                    result.add_genotype_field(sample, "PS", to_string(region_begin(phase.region) + 1));
                    result.add_genotype_field(sample, "PQ", Octopus::to_string(Maths::probability_to_phred<float>(phase.score), 2)); // TODO
                }
            }
        }
        
        return result.build_once();
    }
    
} // namespace Octopus
