// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "vcf_record_factory.hpp"

#include <string>
#include <unordered_set>
#include <set>
#include <unordered_map>
#include <utility>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <sstream>
#include <iostream>

#include <boost/optional.hpp>

#include "concepts/equitable.hpp"
#include "concepts/mappable.hpp"
#include "basics/genomic_region.hpp"
#include "basics/mappable_reference_wrapper.hpp"
#include "core/types/allele.hpp"
#include "core/types/calls/variant_call.hpp"
#include "core/types/calls/call_utils.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/read_stats.hpp"
#include "utils/string_utils.hpp"
#include "utils/maths.hpp"
#include "utils/append.hpp"
#include "exceptions/program_error.hpp"
#include "io/variant/vcf_spec.hpp"
#include "config/octopus_vcf.hpp"

#define _unused(x) ((void)(x))

namespace octopus {

namespace {

constexpr char dummy_base {'#'};
constexpr char deleted_base {'*'};

} // namespace

VcfRecordFactory::VcfRecordFactory(const ReferenceGenome& reference, const ReadMap& reads,
                                   std::vector<SampleName> samples, bool sites_only)
: reference_ {reference}
, reads_ {reads}
, samples_ {std::move(samples)}
, sites_only_ {sites_only}
{}

namespace {

using CallWrapperReference = MappableReferenceWrapper<CallWrapper>;

bool are_in_phase(const Call::GenotypeCall& lhs, const Call::GenotypeCall& rhs)
{
    return lhs.phase && overlaps(lhs.phase->region(), rhs.genotype);
}

} // namespace

class InconsistentCallError : public ProgramError
{
public:
    InconsistentCallError(SampleName sample, Allele first, Allele second)
    : sample_ {std::move(sample)}
    , first_ {std::move(first)}
    , second_ {std::move(second)}
    {}
private:
    SampleName sample_;
    Allele first_, second_;
        
    std::string do_where() const override
    {
        return "VcfRecordFactory::make";
    }
    std::string do_why() const override
    {
        std::ostringstream ss {};
        ss << "In sample " << sample_ << ", alleles " << first_ << " & " << second_ << " were both called";
        return ss.str();
    }
};

void resolve_indel_genotypes(std::vector<CallWrapper>& calls, const std::vector<SampleName>& samples,
                             const ReferenceGenome& reference)
{
    for (auto it = begin(calls); it != end(calls);) {
        if (is_empty(it->mapped_region())) {
            const auto rit = std::find_if_not(make_reverse_iterator(it), make_reverse_iterator(begin(calls)),
                                              [it] (const auto& call) { return are_adjacent(call, *it); });
            // Now everything between rit and it is adjacent to the insertion at it, and will have inserted sequence
            // we want to remove.
            for_each(rit.base(), it, [&samples, &reference, it] (auto& call) {
                for (const auto& sample : samples) {
                    const auto& insertion_genotype = (*it)->get_genotype_call(sample).genotype;
                    auto& sample_genotype = call->get_genotype_call(sample).genotype;
                    std::vector<Allele::NucleotideSequence> resolved_alleles {};
                    resolved_alleles.reserve(sample_genotype.ploidy());
                    transform(std::cbegin(sample_genotype), std::cend(sample_genotype),
                              std::cbegin(insertion_genotype), std::back_inserter(resolved_alleles),
                              [&sample, &reference] (const Allele& allele1, const Allele& allele2) {
                                  if (is_insertion(allele2) && !is_reference(allele1, reference)) {
                                      const auto& old_sequence = allele1.sequence();
                                      if (old_sequence.size() <= sequence_size(allele2)) {
                                          throw InconsistentCallError {sample, allele1, allele2};
                                      }
                                      return Allele::NucleotideSequence {
                                      cbegin(old_sequence), prev(cend(old_sequence), sequence_size(allele2))
                                      };
                                  }
                                  return allele1.sequence();
                              });
                    Genotype<Allele> new_genotype {sample_genotype.ploidy()};
                    for (auto& sequence : resolved_alleles) {
                        new_genotype.emplace(Allele {mapped_region(sample_genotype), move(sequence)});
                    }
                    sample_genotype = std::move(new_genotype);
                }
            });
            auto it2 = std::find_if_not(next(it), end(calls),
                                        [it] (const auto& call) {
                                            return call->mapped_region() == it->mapped_region();
                                        });
            if (it2 == end(calls)) break;
            if (!overlaps(*it, *it2)) {
                it = it2;
                continue;
            }
            auto it3 = find_first_after(next(it2), end(calls), *it);
            // Now everything between it and it2 is an insertion, anything between
            // it2 and it3 is another call which will have inserted sequence we want to remove.
            // Note the genotype calls of all insertions must be the same as they are in the
            // same region
            for_each(it2, it3, [&samples, &reference, it] (auto& call) {
                for (const auto& sample : samples) {
                    const auto& insertion_genotype = (*it)->get_genotype_call(sample).genotype;
                    auto& sample_genotype = call->get_genotype_call(sample).genotype;
                    std::vector<Allele::NucleotideSequence> resolved_alleles {};
                    resolved_alleles.reserve(sample_genotype.ploidy());
                    transform(std::cbegin(sample_genotype), std::cend(sample_genotype),
                              std::cbegin(insertion_genotype), std::back_inserter(resolved_alleles),
                              [&sample, &reference] (const Allele& allele1, const Allele& allele2) {
                                  if (is_insertion(allele2) && !is_reference(allele1, reference)) {
                                      const auto& old_sequence = allele1.sequence();
                                      if (old_sequence.size() <= sequence_size(allele2)) {
                                          throw InconsistentCallError {sample, allele1, allele2};
                                      }
                                      return Allele::NucleotideSequence {
                                      next(cbegin(old_sequence), sequence_size(allele2)), cend(old_sequence)
                                      };
                                  }
                                  return allele1.sequence();
                              });
                    Genotype<Allele> new_genotype {sample_genotype.ploidy()};
                    for (auto& sequence : resolved_alleles) {
                        new_genotype.emplace(Allele {mapped_region(sample_genotype), move(sequence)});
                    }
                    sample_genotype = std::move(new_genotype);
                }
            });
            it = it3;
        } else {
            ++it;
        }
    }
}

bool is_modified_phase_boundary(const CallWrapper& call, const std::vector<SampleName>& samples)
{
    return std::none_of(std::cbegin(samples), std::cend(samples),
                        [&call] (const auto& sample) {
                            const auto& old_phase = call->get_genotype_call(sample).phase;
                            return old_phase && begins_before(call, old_phase->region());
                        });
}

template <typename Iterator>
void resolve_phase(CallWrapper& call, const SampleName& sample,
                   const Iterator first_phase_boundary_itr, const Iterator last_phase_boundary_itr)
{
    const auto& phase = call->get_genotype_call(sample).phase;
    if (phase) {
        auto overlapped = overlap_range(first_phase_boundary_itr, last_phase_boundary_itr, phase->region());
        if (overlapped.empty()) {
            overlapped = overlap_range(first_phase_boundary_itr, last_phase_boundary_itr, expand_lhs(phase->region(), 1));
            if (!overlapped.empty() && begin_distance(overlapped.front(), phase->region()) != 1) {
                overlapped.advance_begin(1);
            }
        }
        if (!overlapped.empty() && overlapped.front() != call) {
            const auto& old_phase = call->get_genotype_call(sample).phase;
            auto new_phase_region = encompassing_region(overlapped.front(), old_phase->region());
            Call::PhaseCall new_phase {std::move(new_phase_region), old_phase->score()};
            call->set_phase(sample, std::move(new_phase));
        }
    }
}

void pad_indels(std::vector<CallWrapper>& calls, const std::vector<SampleName>& samples)
{
    using std::begin; using std::end;
    const auto first_modified_itr = std::stable_partition(begin(calls), end(calls),
                                                          [] (const auto& call) { return !call->parsimonise(dummy_base); });
    if (first_modified_itr != end(calls)) {
        const auto last_call_itr = end(calls);
        const auto first_phase_adjusted_itr = std::partition(first_modified_itr, last_call_itr,
                                                             [&samples] (const auto& call) {
                                                                 return is_modified_phase_boundary(call, samples); });
        if (first_phase_adjusted_itr != last_call_itr) {
            std::sort(first_phase_adjusted_itr, last_call_itr);
            for (auto call_itr = first_phase_adjusted_itr; call_itr != last_call_itr; ++call_itr) {
                auto& call = *call_itr;
                for (const auto& sample : samples) {
                    const auto& old_phase = call->get_genotype_call(sample).phase;
                    if (old_phase) {
                        if (begins_before(call, old_phase->region())) {
                            auto new_phase_region = expand_lhs(old_phase->region(), 1);
                            Call::PhaseCall new_phase {std::move(new_phase_region), old_phase->score()};
                            call->set_phase(sample, std::move(new_phase));
                        } else {
                            resolve_phase(call, sample, first_phase_adjusted_itr, call_itr);
                        }
                    }
                }
            }
            for_each(begin(calls), first_phase_adjusted_itr, [&samples, first_phase_adjusted_itr, last_call_itr] (auto& call) {
                for (const auto& sample : samples) {
                    resolve_phase(call, sample, first_phase_adjusted_itr, last_call_itr);
                }
            });
        }
        std::sort(first_modified_itr, first_phase_adjusted_itr);
        std::inplace_merge(first_modified_itr, first_phase_adjusted_itr, last_call_itr);
        std::inplace_merge(begin(calls), first_modified_itr, last_call_itr);
    }
}

class InconsistentPloidyError : public ProgramError
{
public:
    InconsistentPloidyError(SampleName sample, Genotype<Allele> first, Genotype<Allele> second)
    : sample_ {std::move(sample)}
    , first_ {std::move(first)}
    , second_ {std::move(second)}
    {}
private:
    SampleName sample_;
    Genotype<Allele> first_, second_;
    
    std::string do_where() const override
    {
        return "VcfRecordFactory::make";
    }
    std::string do_why() const override
    {
        std::ostringstream ss {};
        ss << "In sample " << sample_ << ", calls at "
           << first_.mapped_region()
           << " & " << second_.mapped_region()
           << " were called in the same phase set but have different genotype ploidies"
           << " (" << first_.ploidy() << " & "
           << second_.ploidy() << ")";
        return ss.str();
    }
};

std::vector<VcfRecord> VcfRecordFactory::make(std::vector<CallWrapper>&& calls) const
{
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::next;
    using std::prev; using std::for_each; using std::transform; using std::move;
    using std::make_reverse_iterator;
    // TODO: refactor this!!!
    assert(std::is_sorted(std::cbegin(calls), std::cend(calls)));
    resolve_indel_genotypes(calls, samples_, reference_);
    sort_genotype_alleles_by_phase_set(calls, samples_, reference_);
    pad_indels(calls, samples_);
    std::vector<VcfRecord> result {};
    result.reserve(calls.size());
    for (auto call_itr = begin(calls); call_itr != end(calls);) {
        const auto block_begin_itr = adjacent_overlap_find(call_itr, end(calls));
        transform(std::make_move_iterator(call_itr), std::make_move_iterator(block_begin_itr), std::back_inserter(result),
                  [this] (CallWrapper&& call) {
                      call->replace(dummy_base, reference_.fetch_sequence(head_position(call)).front());
                      // We may still have uncalled genotyped alleles here if the called genotype
                      // did not have a high posterior
                      call->replace_uncalled_genotype_alleles(Allele {call->mapped_region(), vcfspec::missingValue}, 'N');
                      return this->make(move(call.call));
                  });
        if (block_begin_itr == end(calls)) break;
        auto block_end_itr = find_next_mutually_exclusive(block_begin_itr, end(calls));
        const auto block_size = std::distance(block_begin_itr, block_end_itr);
        _unused(block_size);
        assert(block_size > 1);
        auto block_head_end_itr = std::find_if_not(next(block_begin_itr), end(calls),
                                                   [block_begin_itr] (const auto& call) {
                                                       return begins_equal(call, *block_begin_itr);
                                                   });
        const auto alt_itr = std::find_if_not(block_begin_itr, block_head_end_itr,
                                              [] (const auto& call) {
                                                  return call->reference().sequence().front() == dummy_base;
                                              });
        boost::optional<decltype(block_head_end_itr)> base {};
        if (alt_itr != block_head_end_itr) base = alt_itr;
        std::deque<CallWrapper> duplicates {};
        for_each(block_begin_itr, block_head_end_itr, [this, base, &duplicates] (auto& call) {
            assert(!call->reference().sequence().empty());
            if (call->reference().sequence().front() == dummy_base) {
                const auto actual_reference_base = reference_.fetch_sequence(head_position(call)).front();
                auto new_sequence = call->reference().sequence();
                new_sequence.front() = actual_reference_base;
                Allele new_allele {mapped_region(call), move(new_sequence)};
                std::unordered_map<Allele, std::set<Allele>> replacements {};
                call->replace(call->reference(), move(new_allele));
                for (const auto& sample : samples_) {
                    auto& genotype = call->get_genotype_call(sample).genotype;
                    for (unsigned i {0}; i < genotype.ploidy(); ++i) {
                        assert(!genotype[i].sequence().empty());
                        if (genotype[i].sequence().front() == dummy_base) {
                            auto new_sequence = genotype[i].sequence();
                            if (base) {
                                const auto& base_genotype = (**base)->get_genotype_call(sample).genotype;
                                if (base_genotype.ploidy() == genotype.ploidy()) {
                                    const auto& base_sequence = base_genotype[i].sequence();
                                    if (!base_sequence.empty()) {
                                        new_sequence.front() = base_sequence.front();
                                    } else {
                                        new_sequence = vcfspec::missingValue;
                                    }
                                } else {
                                    new_sequence.front() = actual_reference_base;
                                }
                            } else {
                                new_sequence.front() = actual_reference_base;
                            }
                            Allele new_allele {mapped_region(call), move(new_sequence)};
                            replacements[genotype[i]].insert(new_allele);
                            genotype[i] = move(new_allele);
                        }
                    }
                }
                for (auto& p : replacements) {
                    std::transform(std::next(std::cbegin(p.second)), std::cend(p.second), std::back_inserter(duplicates),
                                   [&] (const Allele& replacement) {
                                       auto duplicate = clone(call);
                                       duplicate->replace(p.first, replacement);
                                       return duplicate;
                                   });
                    call->replace(p.first, *std::cbegin(p.second));
                }
            }
        });
        if (std::distance(block_begin_itr, block_head_end_itr) > 1) {
            // Here we sort out overlapping records that start at the same position
            auto rit3 = next(make_reverse_iterator(block_head_end_itr));
            const auto rit2 = make_reverse_iterator(block_begin_itr);
            for (; rit3 != rit2; ++rit3) {
                auto& curr_call = *rit3;
                auto& prev_call = *prev(rit3);
                for (const auto& sample : samples_) {
                    auto& curr_genotype = curr_call->get_genotype_call(sample).genotype;
                    auto& prev_genotype = prev_call->get_genotype_call(sample).genotype;
                    for (unsigned i {0}; i < curr_genotype.ploidy(); ++i) {
                        if (prev_genotype.ploidy() > i // can happen if variants not phased but current call padded
                            && (prev_genotype[i].sequence() == vcfspec::deleteMaskAllele ||
                            (curr_genotype[i] != prev_genotype[i]
                             && prev_genotype[i].sequence() == curr_genotype[i].sequence()
                             && sequence_size(curr_genotype[i]) < region_size(curr_genotype)))) {
                            curr_genotype[i] = Allele {mapped_region(curr_call), vcfspec::deleteMaskAllele};
                        } else if (curr_genotype.ploidy() > i
                                && prev_genotype[i] != curr_genotype[i]
                                && sequence_size(curr_genotype[i]) < region_size(curr_genotype) // deletion
                                && sequence_size(prev_genotype[i]) < region_size(prev_genotype) // deletion
                                && prev_genotype[i].sequence().size() > curr_genotype[i].sequence().size()
                                && !prev_call->is_represented(prev_genotype[i])) {
                            prev_genotype[i] = Allele {mapped_region(prev_call), vcfspec::deleteMaskAllele};
                        }
                    }
                }
            }
        }
        std::vector<std::vector<const Call*>> prev_represented {};
        prev_represented.reserve(samples_.size());
        for (const auto& sample : samples_) {
            const auto& genotype = block_begin_itr->call->get_genotype_call(sample).genotype;
            const auto ploidy = genotype.ploidy();
            prev_represented.emplace_back(ploidy, nullptr);
            for (auto itr = block_begin_itr; itr != block_head_end_itr; ++itr) {
                const auto& gt = itr->call->get_genotype_call(sample).genotype;
                if (gt.ploidy() != ploidy) {
                    throw InconsistentPloidyError {sample, genotype, gt};
                }
                for (unsigned i {0}; i < gt.ploidy(); ++i) {
                    if (itr->call->is_represented(gt[i])) {
                        prev_represented.back()[i] = std::addressof(*itr->call);
                    }
                }
            }
        }
        assert(block_begin_itr < block_head_end_itr);
        for (; block_head_end_itr != block_end_itr; ++block_head_end_itr) {
            auto& curr_call = *block_head_end_itr;
            std::unordered_map<Allele, Allele> replacements {};
            assert(!curr_call->reference().sequence().empty());
            if (curr_call->reference().sequence().front() == dummy_base) {
                const auto actual_reference_base = reference_.fetch_sequence(head_position(curr_call)).front();
                auto new_ref_sequence = curr_call->reference().sequence();
                new_ref_sequence.front() = actual_reference_base;
                Allele new_ref_allele {mapped_region(curr_call), move(new_ref_sequence)};
                curr_call->replace(curr_call->reference(), move(new_ref_allele));
                for (unsigned s {0}; s < samples_.size(); ++s) {
                    const auto& sample = samples_[s];
                    auto& genotype_call = curr_call->get_genotype_call(sample);
                    auto& genotype = genotype_call.genotype;
                    for (unsigned i {0}; i < genotype.ploidy(); ++i) {
                        if (genotype[i].sequence().empty()) {
                            Allele::NucleotideSequence new_sequence(region_size(curr_call), deleted_base);
                            genotype[i] = Allele {mapped_region(curr_call), vcfspec::deleteMaskAllele};
                        } else if (genotype[i].sequence().front() == dummy_base) {
                            if (prev_represented[s].size() > i && prev_represented[s][i]
                                && begins_before(*prev_represented[s][i], curr_call)) {
                                const auto& prev_represented_genotype = prev_represented[s][i]->get_genotype_call(sample);
                                if (are_in_phase(genotype_call, prev_represented_genotype)) {
                                    const auto& prev_allele = prev_represented_genotype.genotype[i];
                                    const auto overlap = overlapped_region(prev_allele, curr_call);
                                    if (overlap && prev_allele != prev_represented[s][i]->reference()) {
                                        auto new_sequence = genotype[i].sequence();
                                        const auto overlap_size = static_cast<std::size_t>(region_size(*overlap));
                                        std::fill_n(std::begin(new_sequence), std::min(overlap_size, new_sequence.size()),
                                                    deleted_base);
                                        Allele new_allele {mapped_region(curr_call), move(new_sequence)};
                                        replacements.emplace(genotype[i], new_allele);
                                        genotype[i] = move(new_allele);
                                    } else {
                                        auto new_sequence = genotype[i].sequence();
                                        assert(!genotype[i].sequence().empty());
                                        new_sequence.front() = actual_reference_base;
                                        Allele new_allele {mapped_region(curr_call), move(new_sequence)};
                                        replacements.emplace(genotype[i], new_allele);
                                        genotype[i] = move(new_allele);
                                    }
                                } else {
                                    auto new_sequence = genotype[i].sequence();
                                    assert(!genotype[i].sequence().empty());
                                    new_sequence.front() = actual_reference_base;
                                    Allele new_allele {mapped_region(curr_call), move(new_sequence)};
                                    replacements.emplace(genotype[i], new_allele);
                                    genotype[i] = move(new_allele);
                                }
                            } else {
                                auto new_sequence = genotype[i].sequence();
                                assert(!genotype[i].sequence().empty());
                                new_sequence.front() = actual_reference_base;
                                Allele new_allele {mapped_region(curr_call), move(new_sequence)};
                                replacements.emplace(genotype[i], new_allele);
                                genotype[i] = move(new_allele);
                            }
                        }
                    }
                }
            } else {
                for (unsigned s {0}; s < samples_.size(); ++s) {
                    auto& genotype = curr_call->get_genotype_call(samples_[s]).genotype;
                    for (unsigned i {0}; i < genotype.ploidy(); ++i) {
                        if (genotype[i].sequence().empty()) {
                            Allele::NucleotideSequence new_sequence(region_size(curr_call), deleted_base);
                            genotype[i] = Allele {mapped_region(curr_call), move(new_sequence)};
                        }
                    }
                }
            }
            for (auto& p : replacements) {
                curr_call->replace(p.first, p.second);
            }
            for (unsigned s {0}; s < samples_.size(); ++s) {
                const auto& new_genotype = block_head_end_itr->call->get_genotype_call(samples_[s]).genotype;
                for (unsigned i {0}; i < new_genotype.ploidy(); ++i) {
                    const auto& seq = new_genotype[i].sequence();
                    if (std::find(std::cbegin(seq), std::cend(seq), deleted_base) == std::cend(seq)
                        && block_head_end_itr->call->is_represented(new_genotype[i])) {
                        if (prev_represented[s].size() <= i) {
                            prev_represented[s].resize(i + 1, nullptr);
                        }
                        prev_represented[s][i] = std::addressof(*block_head_end_itr->call);
                    }
                }
            }
        }
        for_each(block_begin_itr, block_end_itr, [] (auto& call) {
            call->replace_uncalled_genotype_alleles(Allele {call->mapped_region(), vcfspec::deleteMaskAllele},
                                                    deleted_base);
        });
        // At this point, all genotypes fields contain canonical bases, '.', or '*', but not '#'.
        std::vector<std::vector<CallWrapper>> segements;
        if (duplicates.empty()) {
            segements = segment_by_begin_copy(std::make_move_iterator(block_begin_itr), std::make_move_iterator(block_end_itr));
        } else {
            std::vector<CallWrapper> section {std::make_move_iterator(block_begin_itr), std::make_move_iterator(block_end_itr)};
            auto itr = utils::append(std::move(duplicates), section);
            std::inplace_merge(std::begin(section), itr, std::end(section));
            segements = segment_by_begin_move(section);
        }
        for (auto&& segment : segements) {
            for (auto&& new_segment : segment_by_end_move(segment)) {
                std::vector<std::unique_ptr<Call>> final_segment {};
                transform(std::make_move_iterator(begin(new_segment)), std::make_move_iterator(end(new_segment)),
                          std::back_inserter(final_segment),
                          [] (auto&& call) -> std::unique_ptr<Call>&& { return move(call.call); });
                result.emplace_back(this->make_segment(move(final_segment)));
            }
        }
        call_itr = block_end_itr;
    }
    return result;
}

// private methods

std::vector<VcfRecord::NucleotideSequence>
extract_all_genotyped_alleles(const Call* call, const std::vector<SampleName>& samples)
{
    using std::begin; using std::end; using std::cbegin; using std::cend;
    std::vector<VcfRecord::NucleotideSequence> result {};
    static constexpr char unique_dummy_base {'~'};
    for (const auto& sample : samples) {
        const auto& called_genotype = call->get_genotype_call(sample).genotype;
        std::transform(cbegin(called_genotype), cend(called_genotype),
                       std::back_inserter(result), [] (const Allele& allele) {
                           auto result = allele.sequence();
                           std::replace(begin(result), end(result), deleted_base, unique_dummy_base);
                           return result;
                       });
    }
    auto itr = std::remove(begin(result), end(result), vcfspec::missingValue);
    result.erase(itr, end(result));
    std::sort(begin(result), end(result));
    itr = std::unique(begin(result), end(result));
    result.erase(itr, end(result));
    for (auto& alt : result) {
        std::replace(begin(alt), end(alt), unique_dummy_base, deleted_base);
    }
    return result;
}

auto extract_genotyped_alt_alleles(const Call* call, const std::vector<SampleName>& samples)
{
    auto result = extract_all_genotyped_alleles(call, samples);
    auto itr = std::find(std::begin(result), std::end(result), call->reference().sequence());
    if (itr != std::end(result)) result.erase(itr);
    return result;
}

auto extract_allele_sequences(const Genotype<Allele>& genotype)
{
    std::vector<VcfRecord::NucleotideSequence> result {};
    result.reserve(genotype.ploidy());
    std::transform(std::cbegin(genotype), std::end(genotype), std::back_inserter(result),
                   [] (const auto& allele) { return allele.sequence(); });
    return result;
}

void set_vcf_genotype(const SampleName& sample, const Call::GenotypeCall& call, VcfRecord::Builder& record,
                      const bool replace_missing_with_non_ref = false)
{
    auto genotyped_alleles = extract_allele_sequences(call.genotype);
    if (replace_missing_with_non_ref) {
        std::replace(std::begin(genotyped_alleles), std::end(genotyped_alleles),
                     std::string {vcfspec::missingValue}, std::string {vcfspec::allele::nonref});
    }
    record.set_genotype(sample, std::move(genotyped_alleles), VcfRecord::Builder::Phasing::phased);
}

bool is_missing(const Allele::NucleotideSequence& sequence) noexcept
{
    return sequence == vcfspec::missingValue;
}
bool is_missing(const Allele& a) noexcept
{
    return is_missing(a.sequence());
}

auto get_allele_counts(const std::vector<VcfRecord::NucleotideSequence>& alt_alleles,
                       const Call& call, const std::vector<SampleName>& samples)
{
    std::vector<unsigned> ac(alt_alleles.size(), 0);
    unsigned an {0};
    for (const auto& sample : samples) {
        for (const auto& allele : call.get_genotype_call(sample).genotype) {
            if (!is_missing(allele)) {
                ++an;
                const auto itr = std::find(std::cbegin(alt_alleles), std::cend(alt_alleles), allele.sequence());
                if (itr != std::cend(alt_alleles)) {
                    ++ac[std::distance(std::cbegin(alt_alleles), itr)];
                }
            }
        }
    }
    return std::make_pair(std::move(ac), an);
}

void set_allele_counts(const Call& call, const std::vector<SampleName>& samples,
                       const std::vector<VcfRecord::NucleotideSequence>& alt_alleles,
                       VcfRecord::Builder& result)
{
    auto p = get_allele_counts(alt_alleles, call, samples);
    result.set_info("AC", utils::to_strings(p.first));
    result.set_info("AN", p.second);
}

auto phred_round(const double phreds)
{
    return phreds < 1 ? maths::round_sf(phreds, 2) : maths::round(phreds, 2);
}

VcfRecord VcfRecordFactory::make(std::unique_ptr<Call> call) const
{
    auto result = VcfRecord::Builder {};
    const auto& region = call->mapped_region();
    bool is_refcall {false};
    auto alts = extract_genotyped_alt_alleles(call.get(), samples_);
    if (alts.empty()) {
        alts.push_back(vcfspec::allele::nonref);
        is_refcall = true;
    } else {
        is_refcall = std::find(std::cbegin(alts), std::cend(alts), vcfspec::allele::nonref) != std::cend(alts);
    }
    result.set_chrom(contig_name(region));
    result.set_pos(mapped_begin(region) + 1);
    result.set_ref(call->reference().sequence());
    set_allele_counts(*call, samples_, alts, result);
    result.set_alt(std::move(alts));
    result.set_qual(std::min(max_qual, phred_round(call->quality().score())));
    const auto call_reads = copy_overlapped(reads_, region);
    result.set_info("NS",  count_samples_with_coverage(call_reads));
    result.set_info("DP",  sum_max_coverages(call_reads));
    result.set_info("MQ",  static_cast<unsigned>(rmq_mapping_quality(call_reads)));
    if (call->model_posterior()) {
        result.set_info("MP",  phred_round(call->model_posterior()->score()));
    }
    if (!sites_only_) {
        std::vector<std::string> format {"GT", "GQ", "DP", "MQ"};
        if (call->all_phased()) {
            utils::append(std::vector<std::string> {"PS", "PQ"}, format);
        }
        const bool add_mp {std::any_of(std::cbegin(samples_), std::cend(samples_), 
                                       [&] (const auto& sample) { return static_cast<bool>(call->model_posterior(sample)); })};
        if (add_mp) format.push_back("MP");
        result.set_format(std::move(format));
        for (const auto& sample : samples_) {
            const auto& genotype_call = call->get_genotype_call(sample);
            static const Phred<double> max_genotype_quality {10'000};
            const auto gq = static_cast<int>(phred_round(std::min(max_genotype_quality, genotype_call.posterior).score()));
            set_vcf_genotype(sample, genotype_call, result, is_refcall);
            result.set_format(sample, "GQ", std::to_string(gq));
            result.set_format(sample, "DP", max_coverage(call_reads.at(sample)));
            result.set_format(sample, "MQ", static_cast<unsigned>(rmq_mapping_quality(call_reads.at(sample))));
            if (call->is_phased(sample)) {
                const auto& phase = *genotype_call.phase;
                auto pq = std::min(100, static_cast<int>(std::round(phase.score().score())));
                result.set_format(sample, "PS", mapped_begin(phase.region()) + 1);
                result.set_format(sample, "PQ", std::to_string(pq));
            }
            if (add_mp) {
                const auto mp = call->model_posterior(sample);
                if (mp) {
                    result.set_format(sample, "MP", phred_round(mp->score()));
                } else {
                    result.set_format_missing(sample, "MP");
                }
            }
        }
    }
    call->decorate(result);
    if (is_refcall) {
        try {
            result.set_blocked_reference();
        } catch (...) {}
    }
    result.collapse_spanning_deletions();
    return result.build_once();
}

namespace {

boost::optional<double> get_model_posterior(const std::vector<std::unique_ptr<Call>>& calls)
{
    std::vector<double> model_posteriors {};
    model_posteriors.reserve(calls.size());
    for (const auto& call : calls) {
        const auto call_model_posterior = call->model_posterior();
        if (call_model_posterior) model_posteriors.push_back(call_model_posterior->score());
    }
    if (model_posteriors.empty()) {
        return boost::none;
    } else {
        return *std::max_element(std::cbegin(model_posteriors), std::cend(model_posteriors));
    }
}

auto get_allele_counts(const std::vector<VcfRecord::NucleotideSequence>& alt_alleles,
                       const std::vector<std::vector<VcfRecord::NucleotideSequence>>& genotypes)
{
    std::vector<unsigned> ac(alt_alleles.size(), 0);
    unsigned an {0};
    for (const auto& genotype : genotypes) {
        for (const auto& allele : genotype) {
            if (!is_missing(allele)) {
                ++an;
                const auto itr = std::find(std::cbegin(alt_alleles), std::cend(alt_alleles), allele);
                if (itr != std::cend(alt_alleles)) {
                    ++ac[std::distance(std::cbegin(alt_alleles), itr)];
                }
            }
        }
    }
    return std::make_pair(std::move(ac), an);
}

void set_allele_counts(const std::vector<VcfRecord::NucleotideSequence>& alt_alleles,
                       const std::vector<std::vector<VcfRecord::NucleotideSequence>>& genotypes,
                       VcfRecord::Builder& result)
{
    auto p = get_allele_counts(alt_alleles, genotypes);
    result.set_info("AC", utils::to_strings(p.first));
    result.set_info("AN", p.second);
}

} // namespace

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
                       [] (const Allele& allele) { return allele.sequence(); });
        std::for_each(std::next(std::cbegin(calls)), std::cend(calls), [&] (const auto& call) {
            const auto& called_genotype = call->get_genotype_call(sample).genotype;
            std::transform(std::cbegin(called_genotype), std::cend(called_genotype),
                           std::cbegin(resolved_sample_genotype), std::begin(resolved_sample_genotype),
                           [&ref] (const Allele& allele, const auto& curr) {
                               const auto& seq = allele.sequence();
                               if (seq.size() < curr.size()
                                   || (!seq.empty() && (seq.front() == '.' || seq.front() == '*' || seq == ref))) {
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
    auto itr = std::remove(std::begin(alt_alleles), std::end(alt_alleles), vcfspec::missingValue);
    itr = std::remove(std::begin(alt_alleles), itr, calls.front()->reference().sequence());
    alt_alleles.erase(itr, std::end(alt_alleles));
    std::sort(std::begin(alt_alleles), std::end(alt_alleles));
    itr = std::unique(std::begin(alt_alleles), std::end(alt_alleles));
    alt_alleles.erase(itr, std::end(alt_alleles));
    bool is_refcall {false};
    if (alt_alleles.empty()) {
        alt_alleles.push_back(vcfspec::allele::nonref);
        is_refcall = true;
    } else {
        is_refcall = std::find(std::cbegin(alt_alleles), std::cend(alt_alleles), vcfspec::allele::nonref) != std::cend(alt_alleles);
    }
    set_allele_counts(alt_alleles, resolved_genotypes, result);
    result.set_alt(std::move(alt_alleles));
    auto q = std::min_element(std::cbegin(calls), std::cend(calls),
                              [] (const auto& lhs, const auto& rhs) { return lhs->quality() < rhs->quality(); });
    result.set_qual(std::min(max_qual, phred_round(q->get()->quality().score())));
    result.set_info("NS", count_samples_with_coverage(reads_, region));
    result.set_info("DP", sum_max_coverages(reads_, region));
    result.set_info("MQ", static_cast<unsigned>(rmq_mapping_quality(reads_, region)));
    const auto site_mp = get_model_posterior(calls);
    if (site_mp) {
        result.set_info("MP", phred_round(*site_mp));
    }
    if (!sites_only_) {
        std::vector<std::string> format {"GT", "GQ", "DP", "MQ"};
        if (calls.front()->all_phased()) {
            utils::append(std::vector<std::string> {"PS", "PQ"}, format);
        }
        const bool add_mp {std::any_of(std::cbegin(samples_), std::cend(samples_), 
                                       [&] (const auto& sample) { return static_cast<bool>(calls.front()->model_posterior(sample)); })};
        if (add_mp) format.push_back("MP");
        result.set_format(std::move(format));
        auto sample_itr = std::begin(resolved_genotypes);
        for (const auto& sample : samples_) {
            const auto posterior = calls.front()->get_genotype_call(sample).posterior;
            static const Phred<double> max_genotype_quality {10'000};
            const auto gq = static_cast<int>(phred_round(std::min(max_genotype_quality, posterior).score()));
            auto& genotype_call = *sample_itr++;
            if (is_refcall) {
                std::replace(std::begin(genotype_call), std::end(genotype_call),
                             std::string {vcfspec::missingValue}, std::string {vcfspec::allele::nonref});
            }
            result.set_genotype(sample, genotype_call, VcfRecord::Builder::Phasing::phased);
            result.set_format(sample, "GQ", std::to_string(gq));
            result.set_format(sample, "DP", max_coverage(reads_.at(sample), region));
            result.set_format(sample, "MQ", static_cast<unsigned>(rmq_mapping_quality(reads_.at(sample), region)));
            if (calls.front()->is_phased(sample)) {
                const auto phase = *calls.front()->get_genotype_call(sample).phase;
                auto pq = std::min(100, static_cast<int>(std::round(phase.score().score())));
                result.set_format(sample, "PS", mapped_begin(phase.region()) + 1);
                result.set_format(sample, "PQ", std::to_string(pq));
            }
            if (add_mp) {
                const auto mp = calls.front()->model_posterior(sample);
                if (mp) {
                    result.set_format(sample, "MP", phred_round(mp->score()));
                } else {
                    result.set_format_missing(sample, "MP");
                }
            }
        }
    }
    for (const auto& call : calls) {
        call->decorate(result);
    }
    if (is_refcall) {
        try {
            result.set_blocked_reference();
        } catch (...) {}
    }
    result.collapse_spanning_deletions();
    return result.build_once();
}
    
} // namespace octopus
