// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

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

#include "concepts/equitable.hpp"
#include "basics/genomic_region.hpp"
#include "concepts/mappable.hpp"
#include "utils/mappable_algorithms.hpp"
#include "core/types/allele.hpp"
#include "utils/read_stats.hpp"
#include "utils/string_utils.hpp"
#include "utils/maths.hpp"
#include "variant_call.hpp"

#include "timers.hpp"

namespace octopus {

VcfRecordFactory::VcfRecordFactory(const ReferenceGenome& reference, const ReadMap& reads,
                                   std::vector<SampleName> samples, bool sites_only)
: reference_ {reference}
, reads_ {reads}
, samples_ {std::move(samples)}
, sites_only_ {sites_only}
{}

namespace {

struct CallWrapper : public Equitable<CallWrapper>, public Mappable<CallWrapper>
{
    CallWrapper(std::unique_ptr<Call> call) : call {std::move(call) } {}
    
    CallWrapper(const CallWrapper&)            = delete;
    CallWrapper& operator=(const CallWrapper&) = delete;
    CallWrapper(CallWrapper&&)                 = default;
    CallWrapper& operator=(CallWrapper&&)      = default;
    
    operator const std::unique_ptr<Call>&() const noexcept { return call; }
    operator std::unique_ptr<Call>&() noexcept { return call; }
    std::unique_ptr<Call>::pointer operator->() const noexcept { return call.get(); };
    const GenomicRegion& mapped_region() const noexcept { return call->mapped_region(); }
    
    std::unique_ptr<Call> call;
};

bool operator==(const CallWrapper& lhs, const CallWrapper& rhs)
{
    return *lhs.call == *rhs.call;
}

template <typename T>
auto wrap(std::vector<std::unique_ptr<T>>&& calls)
{
    return std::vector<CallWrapper> {
        std::make_move_iterator(std::begin(calls)),
        std::make_move_iterator(std::end(calls))
    };
}

struct CallWrapperReference : public Mappable<CallWrapperReference>
{
    std::reference_wrapper<const CallWrapper> call;
    
    template <typename T> CallWrapperReference(T&& call) : call {std::forward<T>(call)} {}
    const GenomicRegion& mapped_region() const noexcept { return call.get()->mapped_region(); }
    decltype(auto) get() const noexcept { return call.get(); }
};

} // namespace

bool are_in_phase(const Call::GenotypeCall& lhs, const Call::GenotypeCall& rhs)
{
    return lhs.phase && overlaps(lhs.phase->region(), rhs.genotype);
}

std::vector<VcfRecord> VcfRecordFactory::make(std::vector<std::unique_ptr<Call>>&& calls) const
{
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::next;
    using std::prev; using std::for_each; using std::transform; using std::move;
    
    // TODO: refactor this!!!
    
    auto wrapped_calls = wrap(move(calls));
    
    calls.clear();
    calls.shrink_to_fit();
    
    for (auto it = begin(wrapped_calls); it != end(wrapped_calls);) {
        if (is_empty(it->mapped_region())) {
            auto it2 = std::find_if_not(next(it), end(wrapped_calls),
                                        [it] (const auto& call) {
                                            return call->mapped_region() == it->mapped_region();
                                        });
            if (it2 == end(wrapped_calls)) break;
            if (!overlaps(*it, *it2)) {
                it = it2;
                continue;
            }
            auto it3 = find_first_after(next(it2), end(wrapped_calls), *it);
            
            // Now everything between it and it2 is an insertion, anything between
            // it2 and it3 is another call which will have inserted sequence we want to remove.
            // Note the genotype calls of all insertions must be the same as they are in the
            // same region
            for_each(it2, it3, [this, it] (auto& call) {
                for (const auto& sample : samples_) {
                    const auto& insertion_genotype = (*it)->get_genotype_call(sample).genotype;
                    auto& sample_genotype = call->get_genotype_call(sample).genotype;
                    std::vector<Allele::NucleotideSequence> resolved_alleles {};
                    resolved_alleles.reserve(sample_genotype.ploidy());
                    
                    transform(cbegin(sample_genotype), cend(sample_genotype),
                              cbegin(insertion_genotype), std::back_inserter(resolved_alleles),
                              [] (const Allele& allele1, const Allele& allele2) {
                                  if (is_insertion(allele2)) {
                                      const auto& old_sequence = allele1.sequence();
                                      return Allele::NucleotideSequence {
                                          next(cbegin(old_sequence), sequence_size(allele2)),
                                               cend(old_sequence)
                                      };
                                  }
                                  return allele1.sequence();
                              });
                    
                    Genotype<Allele> new_genotype {sample_genotype.ploidy()};
                    for (auto& sequence : resolved_alleles) {
                        new_genotype.emplace(Allele {mapped_region(sample_genotype), move(sequence)});
                    }
                    sample_genotype = move(new_genotype);
                }
            });
            it = it3;
        } else {
            ++it;
        }
    }
    
    const auto first_modified = std::stable_partition(begin(wrapped_calls), end(wrapped_calls),
                                                      [] (const auto& call) {
                                                          return !call->parsimonise('#');
                                                      });
    if (first_modified != end(wrapped_calls)) {
        const auto last = end(wrapped_calls);
        const auto first_phase_adjusted = std::partition(first_modified, last,
                                        [this] (const auto& call) {
                                            return std::none_of(begin(samples_), cend(samples_),
                                                                [&call] (const auto& sample) {
                                                                    const auto& old_phase = call->get_genotype_call(sample).phase;
                                                                    return old_phase && begins_before(mapped_region(call), old_phase->region());
                                                                });
                                        });
        if (first_phase_adjusted != last) {
            std::sort(first_phase_adjusted, last);
            for_each(first_phase_adjusted, last,
                     [this] (auto& call) {
                         for (const auto& sample : samples_) {
                             const auto& old_phase = call->get_genotype_call(sample).phase;
                             if (old_phase && begins_before(mapped_region(call), old_phase->region())) {
                                 auto new_phase_region = expand_lhs(old_phase->region(), 1);
                                 Call::PhaseCall new_phase {move(new_phase_region), old_phase->score()};
                                 call->set_phase(sample, move(new_phase));
                             }
                         }
                     });
            for_each(begin(wrapped_calls), first_phase_adjusted,
                     [this, first_phase_adjusted, last] (auto& call) {
                         for (const auto& sample : samples_) {
                             const auto& phase = call->get_genotype_call(sample).phase;
                             if (phase) {
                                 auto overlapped = overlap_range(first_phase_adjusted, last, phase->region());
                                 if (overlapped.empty()) {
                                     overlapped = overlap_range(first_phase_adjusted, last, expand_lhs(phase->region(), 1));
                                     if (!overlapped.empty()) {
                                         if (begin_distance(overlapped.front(), phase->region()) != 1) {
                                             overlapped.advance_begin(1);
                                         }
                                     }
                                 }
                                 if (!overlapped.empty() && overlapped.front() != call) {
                                     const auto& old_phase = call->get_genotype_call(sample).phase;
                                     auto new_phase_region = encompassing_region(overlapped.front(), old_phase->region());
                                     Call::PhaseCall new_phase {move(new_phase_region), old_phase->score()};
                                     call->set_phase(sample, move(new_phase));
                                 }
                             }
                         }
                     });
        }
        std::sort(first_modified, first_phase_adjusted);
        std::inplace_merge(first_modified, first_phase_adjusted, last);
        std::inplace_merge(begin(wrapped_calls), first_modified, last);
    }
    
    std::vector<VcfRecord> result {};
    result.reserve(calls.size());
    
    for (auto it = begin(wrapped_calls); it != end(wrapped_calls);) {
        const auto it2 = adjacent_overlap_find(it, end(wrapped_calls));
        transform(std::make_move_iterator(it), std::make_move_iterator(it2),
                  std::back_inserter(result),
                  [this] (CallWrapper&& call) {
                      call->replace('#', reference_.fetch_sequence(head_position(call->mapped_region())).front());
                      // We may still have uncalled genotyped alleles here if the called genotype
                      // did not have a high posterior
                      call->replace_uncalled_genotype_alleles(Allele {call->mapped_region(), "."}, 'N');
                      return this->make(move(call.call));
                  });
        if (it2 == end(wrapped_calls)) break;
        auto it3 = std::find_if_not(next(it2), end(wrapped_calls),
                                    [it2] (const auto& call) {
                                        return begins_equal(call, *it2);
                                    });
        boost::optional<decltype(it3)> base;
        if ((*it2)->reference().sequence().front() != '#') {
            base = it2;
        }
        for_each(base ? next(it2) : it2, it3, [this, base] (auto& call) {
            if (call->reference().sequence().front() == '#') {
                const auto actual_reference_base = reference_.fetch_sequence(head_position(call)).front();
                auto new_sequence = call->reference().sequence();
                new_sequence.front() = actual_reference_base;
                Allele new_allele {mapped_region(call), move(new_sequence)};
                std::unordered_map<Allele, Allele> replacements {};
                
                call->replace(call->reference(), move(new_allele));
                for (const auto& sample : samples_) {
                    auto& genotype_call = call->get_genotype_call(sample);
                    auto& old_genotype = genotype_call.genotype;
                    const auto ploidy = old_genotype.ploidy();
                    Genotype<Allele> new_genotype {ploidy};
                    
                    for (unsigned i {0}; i < ploidy; ++i) {
                        if (old_genotype[i].sequence().front() == '#') {
                            auto new_sequence = old_genotype[i].sequence();
                            if (base) {
                                const auto& base_sequence = (**base)->get_genotype_call(sample).genotype[i].sequence();
                                new_sequence.front() = base_sequence.front();
                            } else {
                                new_sequence.front() = actual_reference_base;
                            }
                            Allele new_allele {mapped_region(call), move(new_sequence)};
                            replacements.emplace(old_genotype[i], new_allele);
                            new_genotype.emplace(move(new_allele));
                        } else {
                            new_genotype.emplace(old_genotype[i]);
                        }
                    }
                    old_genotype = move(new_genotype);
                }
                for (auto& p : replacements) {
                    call->replace(p.first, p.second);
                }
            }
        });
        if (std::distance(it2, it3) > 1) {
            auto rit3 = next(std::make_reverse_iterator(it3));
            const auto rit2 = std::make_reverse_iterator(it2);
            for (; rit3 != rit2; ++rit3) {
                auto& curr_call = *rit3;
                for (const auto& sample : samples_) {
                    auto& genotype_call = curr_call->get_genotype_call(sample);
                    auto& old_genotype = genotype_call.genotype;
                    const auto& prev_call = *prev(rit3);
                    const auto& prev_genotype_call = prev_call->get_genotype_call(sample);
                    const auto& prev_genotype = prev_genotype_call.genotype;
                    const auto ploidy = old_genotype.ploidy();
                    Genotype<Allele> new_genotype {ploidy};
                    
                    for (unsigned i {0}; i < ploidy; ++i) {
                        if (prev_genotype[i].sequence() == "*" ||
                            (prev_genotype[i].sequence() == old_genotype[i].sequence()
                             && sequence_size(old_genotype[i]) < region_size(old_genotype))) {
                            Allele::NucleotideSequence new_sequence(1, '*');
                            Allele new_allele {mapped_region(curr_call), move(new_sequence)};
                            new_genotype.emplace(move(new_allele));
                        } else {
                            new_genotype.emplace(old_genotype[i]);
                        }
                    }
                    old_genotype = move(new_genotype);
                }
            }
        }
        if (it3 != end(wrapped_calls) && overlaps(*prev(it3), *it3)) {
            it = find_first_after(next(it3), end(wrapped_calls), *prev(it3));
            while (it != end(wrapped_calls)) {
                if (overlaps(*it, *prev(it))) {
                    it = find_first_after(next(it), end(wrapped_calls), *it);
                } else {
                    break;
                }
            }
        } else {
            it = it3;
        }
        for (; it3 != it; ++it3) {
            auto& curr_call = *it3;
            std::unordered_map<Allele, Allele> replacements {};
            if (curr_call->reference().sequence().front() == '#') {
                const auto actual_reference_base = reference_.fetch_sequence(head_position(curr_call)).front();
                auto new_sequence = curr_call->reference().sequence();
                new_sequence.front() = actual_reference_base;
                Allele new_allele {mapped_region(curr_call), move(new_sequence)};
                curr_call->replace(curr_call->reference(), move(new_allele));
                const auto& prev_call = *prev(it3);
                
                if (overlaps(curr_call, prev_call)) {
                    for (const auto& sample : samples_) {
                        auto& genotype_call = curr_call->get_genotype_call(sample);
                        auto& old_genotype = genotype_call.genotype;
                        const auto& prev_genotype_call = prev_call->get_genotype_call(sample);
                        const auto& prev_genotype = prev_genotype_call.genotype;
                        const auto ploidy = old_genotype.ploidy();
                        assert(ploidy == prev_genotype.ploidy());
                        Genotype<Allele> new_genotype {ploidy};
                        const auto& call_region = mapped_region(curr_call);
                        
                        if (are_in_phase(prev_genotype_call, genotype_call)) {
                            for (unsigned i {0}; i < ploidy; ++i) {
                                if (old_genotype[i].sequence().empty()) {
                                    Allele::NucleotideSequence new_sequence(region_size(curr_call), '*');
                                    Allele new_allele {call_region, move(new_sequence)};
                                    new_genotype.emplace(move(new_allele));
                                } else if (old_genotype[i].sequence().front() == '#') {
                                    if (splice(prev_genotype[i], head_position(curr_call)).sequence().front() == actual_reference_base) {
                                        auto new_sequence = old_genotype[i].sequence();
                                        new_sequence.front() = actual_reference_base;
                                        Allele new_allele {call_region, move(new_sequence)};
                                        replacements.emplace(old_genotype[i], new_allele);
                                        new_genotype.emplace(move(new_allele));
                                    } else {
                                        auto new_sequence = old_genotype[i].sequence();
                                        new_sequence.front() = '*';
                                        Allele new_allele {call_region, move(new_sequence)};
                                        new_genotype.emplace(move(new_allele));
                                    }
                                } else {
                                    new_genotype.emplace(old_genotype[i]);
                                }
                            }
                        } else {
                            for (unsigned i {0}; i < ploidy; ++i) {
                                if (old_genotype[i].sequence().empty()) {
                                    Allele::NucleotideSequence new_sequence(region_size(curr_call), '*');
                                    Allele new_allele {call_region, move(new_sequence)};
                                    new_genotype.emplace(move(new_allele));
                                } else if (old_genotype[i].sequence().front() == '#') {
                                    auto new_sequence = old_genotype[i].sequence();
                                    new_sequence.front() = '.';
                                    Allele new_allele {call_region, move(new_sequence)};
                                    new_genotype.emplace(move(new_allele));
                                } else {
                                    new_genotype.emplace(old_genotype[i]);
                                }
                            }
                        }
                        old_genotype = move(new_genotype);
                    }
                } else {
                    for (const auto& sample : samples_) {
                        auto& genotype_call = curr_call->get_genotype_call(sample);
                        auto& old_genotype = genotype_call.genotype;
                        const auto ploidy = old_genotype.ploidy();
                        Genotype<Allele> new_genotype {ploidy};
                        for (unsigned i {0}; i < ploidy; ++i) {
                            if (old_genotype[i].sequence().front() == '#') {
                                auto new_sequence = old_genotype[i].sequence();
                                new_sequence.front() = actual_reference_base;
                                Allele new_allele {mapped_region(curr_call), move(new_sequence)};
                                replacements.emplace(old_genotype[i], new_allele);
                                new_genotype.emplace(move(new_allele));
                            } else {
                                new_genotype.emplace(old_genotype[i]);
                            }
                        }
                        old_genotype = move(new_genotype);
                    }
                }
            } else {
                for (const auto& sample : samples_) {
                    auto& genotype_call = curr_call->get_genotype_call(sample);
                    auto& old_genotype = genotype_call.genotype;
                    const auto ploidy = old_genotype.ploidy();
                    Genotype<Allele> new_genotype {ploidy};
                    for (unsigned i {0}; i < ploidy; ++i) {
                        if (old_genotype[i].sequence().empty()) {
                            Allele::NucleotideSequence new_sequence(region_size(curr_call), '*');
                            Allele new_allele {mapped_region(curr_call), move(new_sequence)};
                            new_genotype.emplace(move(new_allele));
                        } else {
                            new_genotype.emplace(old_genotype[i]);
                        }
                    }
                    old_genotype = move(new_genotype);
                }
            }
            for (auto& p : replacements) {
                curr_call->replace(p.first, p.second);
            }
        }
        
        for_each(it2, it, [] (auto& call) {
            call->replace_uncalled_genotype_alleles(Allele {call->mapped_region(), "."}, '*');
        });
        
        // At this point, all genotypes fields contain canonical bases, '.', or '*', but not '#'.
        auto segements = segment_by_begin_copy(std::make_move_iterator(it2), std::make_move_iterator(it));
        for (auto&& segment : segements) {
            for (auto&& new_segment : segment_by_end_move(segment)) {
                std::vector<std::unique_ptr<Call>> final_segment {};
                transform(std::make_move_iterator(begin(new_segment)),
                          std::make_move_iterator(end(new_segment)),
                          std::back_inserter(final_segment),
                          [] (auto&& call) -> std::unique_ptr<Call>&& {
                              return move(call.call);
                          });
                result.emplace_back(this->make_segment(move(final_segment)));
            }
        }
    }
    
    return result;
}

// private methods

std::vector<VcfRecord::NucleotideSequence>
extract_all_genotyped_alleles(const Call* call, const std::vector<SampleName>& samples)
{
    using std::begin; using std::end; using std::cbegin; using std::cend;
    
    std::vector<VcfRecord::NucleotideSequence> result {};
    
    for (const auto& sample : samples) {
        const auto& called_genotype = call->get_genotype_call(sample).genotype;
        std::transform(cbegin(called_genotype), cend(called_genotype),
                       std::back_inserter(result), [] (const Allele& allele) {
                           auto result = allele.sequence();
                           std::replace(begin(result), end(result), '*', '~');
                           return result;
                       });
    }
    
    const VcfRecord::NucleotideSequence missing_allele(1, '.');
    auto it = std::remove(begin(result), end(result), missing_allele);
    result.erase(it, end(result));
    std::sort(begin(result), end(result));
    it = std::unique(begin(result), end(result));
    result.erase(it, end(result));
    
    for (auto& alt : result) {
        std::replace(begin(alt), end(alt), '~', '*');
    }
    
    return result;
}

void set_alt_alleles(const Call* call, VcfRecord::Builder& record,
                     const std::vector<SampleName>& samples)
{
    auto alts = extract_all_genotyped_alleles(call, samples);
    auto it = std::find(std::begin(alts), std::end(alts), call->reference().sequence());
    if (it != std::end(alts)) alts.erase(it);
    assert(std::find(std::cbegin(alts), std::cend(alts), "#") == std::cend(alts));
    assert(std::find(std::cbegin(alts), std::cend(alts), ".") == std::cend(alts));
    assert(std::find(std::cbegin(alts), std::cend(alts), "") == std::cend(alts));
    record.set_alt(std::move(alts));
}

void set_vcf_genotype(const SampleName& sample, const Call::GenotypeCall& genotype_call,
                      VcfRecord::Builder& record)
{
    std::vector<VcfRecord::NucleotideSequence> result {};
    result.reserve(genotype_call.genotype.ploidy());
    for (const auto& allele : genotype_call.genotype) {
        result.push_back(allele.sequence());
    }
    record.set_genotype(sample, result, VcfRecord::Builder::Phasing::phased);
}

VcfRecord VcfRecordFactory::make(std::unique_ptr<Call> call) const
{
    auto result = VcfRecord::Builder {};
    const auto& region = call->mapped_region();
    
    result.set_chrom(contig_name(region));
    result.set_pos(mapped_begin(region) + 1);
    result.set_ref(call->reference().sequence());
    set_alt_alleles(call.get(), result, samples_);
    result.set_qual(std::min(5000.0, maths::round(call->quality().score(), 2)));
    const auto call_reads = copy_overlapped(reads_, region);
    result.set_info("NS",  count_samples_with_coverage(call_reads));
    result.set_info("DP",  sum_max_coverages(call_reads));
    result.set_info("SB",  utils::to_string(strand_bias(call_reads), 2));
    result.set_info("BQ",  static_cast<unsigned>(rmq_base_quality(call_reads)));
    result.set_info("MQ",  static_cast<unsigned>(rmq_mapping_quality(call_reads)));
    result.set_info("MQ0", count_mapq_zero(call_reads));
    
    if (call->model_posterior()) {
        result.set_info("MP",  maths::round(call->model_posterior()->score(), 2));
    }
    if (!sites_only_) {
        if (call->all_phased()) {
            result.set_format({"GT", "GQ", "DP", "BQ", "MQ", "PS", "PQ"});
        } else {
            result.set_format({"GT", "GQ", "DP", "BQ", "MQ"});
        }
        
        for (const auto& sample : samples_) {
            const auto& genotype_call = call->get_genotype_call(sample);
            auto gq = std::min(99, static_cast<int>(std::round(genotype_call.posterior.score())));
            
            set_vcf_genotype(sample, genotype_call, result);
            result.set_format(sample, "GQ", std::to_string(gq));
            result.set_format(sample, "DP", max_coverage(call_reads.at(sample)));
            result.set_format(sample, "BQ", static_cast<unsigned>(rmq_base_quality(call_reads.at(sample))));
            result.set_format(sample, "MQ", static_cast<unsigned>(rmq_mapping_quality(call_reads.at(sample))));
            if (call->is_phased(sample)) {
                const auto& phase = *genotype_call.phase;
                auto pq = std::min(99, static_cast<int>(std::round(phase.score().score())));
                result.set_format(sample, "PS", mapped_begin(phase.region()) + 1);
                result.set_format(sample, "PQ", std::to_string(pq));
            }
        }
    }
    
    call->decorate(result);
    
    return result.build_once();
}

VcfRecord VcfRecordFactory::make_segment(std::vector<std::unique_ptr<Call>>&& calls) const
{
    assert(!calls.empty());
    
    if (calls.size() == 1) {
        return make(std::move(calls.front()));
    }
    
    auto result = VcfRecord::Builder {};
    const auto& region = calls.front()->mapped_region();
    const auto& ref = calls.front()->reference().sequence();
    
    result.set_chrom(contig_name(region));
    result.set_pos(mapped_begin(region) + 1);
    result.set_ref(ref);
    
    std::vector<std::vector<VcfRecord::NucleotideSequence>> resolved_genotypes {};
    resolved_genotypes.reserve(samples_.size());
    
    for (const auto& sample : samples_) {
        const auto ploidy = calls.front()->get_genotype_call(sample).genotype.ploidy();
        std::vector<VcfRecord::NucleotideSequence> resolved_sample_genotype(ploidy);
        const auto& first_called_genotype = calls.front()->get_genotype_call(sample).genotype;
        
        std::transform(std::cbegin(first_called_genotype), std::cend(first_called_genotype),
                       std::begin(resolved_sample_genotype),
                       [] (const Allele& allele) {
                           return allele.sequence();
                       });
        std::for_each(std::next(std::cbegin(calls)), std::cend(calls), [&] (const auto& call) {
            const auto& called_genotype = call->get_genotype_call(sample).genotype;
            std::transform(std::cbegin(called_genotype), std::cend(called_genotype),
                           std::cbegin(resolved_sample_genotype),
                           std::begin(resolved_sample_genotype),
                           [&ref] (const Allele& allele, const auto& curr) {
                               const auto& seq = allele.sequence();
                               if (seq.front() == '.' || seq.front() == '*' || seq == ref) {
                                   return curr;
                               }
                               return seq;
                           });
        });
        
        resolved_genotypes.push_back(std::move(resolved_sample_genotype));
    }
    
    std::vector<VcfRecord::NucleotideSequence> alt_alleles {};
    for (const auto& genotype : resolved_genotypes) {
        alt_alleles.insert(std::end(alt_alleles), std::cbegin(genotype), std::cend(genotype));
    }
    const VcfRecord::NucleotideSequence missing_allele(1, '.');
    auto it = std::remove(std::begin(alt_alleles), std::end(alt_alleles), missing_allele);
    it = std::remove(std::begin(alt_alleles), it, calls.front()->reference().sequence());
    alt_alleles.erase(it, std::end(alt_alleles));
    std::sort(std::begin(alt_alleles), std::end(alt_alleles));
    it = std::unique(std::begin(alt_alleles), std::end(alt_alleles));
    alt_alleles.erase(it, std::end(alt_alleles));
    result.set_alt(std::move(alt_alleles));
    auto q = std::max_element(std::cbegin(calls), std::cend(calls),
                              [] (const auto& lhs, const auto& rhs) {
                                  return lhs->quality() < rhs->quality();
                              });
    result.set_qual(std::min(5000.0, maths::round(q->get()->quality().score(), 2)));
    result.set_info("NS",  count_samples_with_coverage(reads_, region));
    result.set_info("DP",  sum_max_coverages(reads_, region));
    result.set_info("SB",  utils::to_string(strand_bias(reads_, region), 2));
    result.set_info("BQ",  static_cast<unsigned>(rmq_base_quality(reads_, region)));
    result.set_info("MQ",  static_cast<unsigned>(rmq_mapping_quality(reads_, region)));
    result.set_info("MQ0", count_mapq_zero(reads_, region));
    
    boost::optional<double> mp {};
    for (const auto& call : calls) {
        const auto call_mp = call->model_posterior();
        if (call_mp) {
            if (mp) {
                if (*mp < call_mp->score()) mp = call_mp->score();
            } else {
                mp = call_mp->score();
            }
        }
    }
    if (mp) {
        result.set_info("MP", maths::round(*mp, 2));
    }
    if (!sites_only_) {
        if (calls.front()->all_phased()) {
            result.set_format({"GT", "GQ", "DP", "BQ", "MQ", "PS", "PQ"});
        } else {
            result.set_format({"GT", "GQ", "DP", "BQ", "MQ"});
        }
        
        auto sample_itr = std::begin(resolved_genotypes);
        for (const auto& sample : samples_) {
            const auto posterior = calls.front()->get_genotype_call(sample).posterior;
            auto gq = std::min(99, static_cast<int>(std::round(posterior.score())));
            
            result.set_genotype(sample, *sample_itr++, VcfRecord::Builder::Phasing::phased);
            result.set_format(sample, "GQ", std::to_string(gq));
            result.set_format(sample, "DP", max_coverage(reads_.at(sample), region));
            result.set_format(sample, "BQ", static_cast<unsigned>(rmq_base_quality(reads_.at(sample), region)));
            result.set_format(sample, "MQ", static_cast<unsigned>(rmq_mapping_quality(reads_.at(sample), region)));
            if (calls.front()->is_phased(sample)) {
                const auto phase = *calls.front()->get_genotype_call(sample).phase;
                auto pq = std::min(99, static_cast<int>(std::round(phase.score().score())));
                result.set_format(sample, "PS", mapped_begin(phase.region()) + 1);
                result.set_format(sample, "PQ", std::to_string(pq));
            }
        }
    }
    
    calls.front()->decorate(result);
    
    return result.build_once();
}
    
} // namespace octopus
