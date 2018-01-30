// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "genotype_reader.hpp"

#include <utility>
#include <functional>
#include <iterator>
#include <algorithm>
#include <initializer_list>
#include <cassert>

#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>

#include "basics/contig_region.hpp"
#include "basics/genomic_region.hpp"
#include "concepts/mappable.hpp"
#include "utils/mappable_algorithms.hpp"
#include "core/types/allele.hpp"
#include "core/types/variant.hpp"
#include "io/variant/vcf_header.hpp"
#include "io/variant/vcf_spec.hpp"
#include "io/reference/reference_genome.hpp"

namespace octopus {

namespace {

bool is_missing(const VcfRecord::NucleotideSequence& allele) noexcept
{
    return allele == vcfspec::missingValue;
}

bool is_delete_masked(const VcfRecord::NucleotideSequence& allele) noexcept
{
    const std::string deleted_sequence {vcfspec::deletedBase};
    return allele == deleted_sequence;
}

bool is_partially_delete_masked(const VcfRecord::NucleotideSequence& allele) noexcept
{
    return allele.size() > 1 && allele.front() == vcfspec::deletedBase && allele.back() != vcfspec::deletedBase;
}

void remove_missing_alleles(std::vector<VcfRecord::NucleotideSequence>& genotype)
{
    genotype.erase(std::remove_if(std::begin(genotype), std::end(genotype), is_missing), std::end(genotype));
}

void remove_deleted_alleles(std::vector<VcfRecord::NucleotideSequence>& genotype)
{
    genotype.erase(std::remove_if(std::begin(genotype), std::end(genotype), is_delete_masked), std::end(genotype));
}

bool is_complex(const VcfRecord::NucleotideSequence& ref, const VcfRecord::NucleotideSequence& alt) noexcept
{
    return !ref.empty() && !alt.empty() && ref.size() != alt.size() && ref.front() != alt.front();
}

bool is_ref_pad_size_known(const VcfRecord::NucleotideSequence& allele, const VcfRecord& call) noexcept
{
    return allele != call.ref() && !is_complex(call.ref(), allele);
}

auto num_matching_lhs_bases(const VcfRecord::NucleotideSequence& lhs, const VcfRecord::NucleotideSequence& rhs) noexcept
{
    auto p = std::mismatch(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs));
    return static_cast<int>(std::distance(std::cbegin(lhs), p.first));
}

auto calculate_ref_pad_size(const VcfRecord& call, const VcfRecord::NucleotideSequence& allele) noexcept
{
    if (allele == "*") {
        return 1;
    } else if (is_partially_delete_masked(allele)) {
        auto first_base_itr = std::find_if_not(std::cbegin(allele), std::cend(allele),
                                               [] (auto b) { return b == vcfspec::deletedBase; });
        return static_cast<int>(std::distance(std::cbegin(allele), first_base_itr));
    } else {
        return num_matching_lhs_bases(call.ref(), allele);
    }
}

bool has_simple_indel(const VcfRecord& call) noexcept
{
    return std::any_of(std::cbegin(call.alt()), std::cend(call.alt()),
                       [&] (const auto& allele) {
                           return allele.size() != call.ref().size() && (allele.size() == 1 || call.ref().size() == 1);
                       });
}

boost::optional<ContigAllele> make_allele(const VcfRecord& call, VcfRecord::NucleotideSequence allele_sequence, const int ref_pad)
{
    if (is_missing(allele_sequence)) {
        return boost::none;
    }
    auto region = contig_region(call);
    if (is_delete_masked(allele_sequence)) {
        if (call.alt().size() > 1) {
            return boost::none;
        } else {
            allele_sequence.clear();
            region = expand_lhs(region, -1);
        }
    } else if (is_partially_delete_masked(allele_sequence)) {
        auto first_base_itr = std::find_if_not(std::cbegin(allele_sequence), std::cend(allele_sequence),
                                               [] (auto b) { return b == vcfspec::deletedBase; });
        auto delete_mask_len = std::distance(std::cbegin(allele_sequence), first_base_itr);
        allele_sequence.erase(std::cbegin(allele_sequence), first_base_itr);
        region = expand_lhs(region, -delete_mask_len);
    } else if (ref_pad > 0) {
        assert(static_cast<std::size_t>(ref_pad) <= allele_sequence.size());
        allele_sequence.erase(std::cbegin(allele_sequence), std::next(std::cbegin(allele_sequence), ref_pad));
        region = expand_lhs(region, -ref_pad);
    }
    return ContigAllele {region, std::move(allele_sequence)};
}

auto extract_genotype(const VcfRecord& call, const SampleName& sample)
{
    auto genotype = get_genotype(call, sample);
    boost::optional<int> min_ref_pad {};
    std::vector<std::size_t> unknown_pad_indices {};
    const auto ploidy = genotype.size();
    std::vector<boost::optional<ContigAllele>> result(ploidy, boost::none);
    for (std::size_t i {0}; i < ploidy; ++i) {
        auto& allele = genotype[i];
        if (is_ref_pad_size_known(allele, call)) {
            const auto allele_pad = num_matching_lhs_bases(call.ref(), allele);
            if (min_ref_pad) {
                min_ref_pad = std::min(*min_ref_pad, allele_pad);
            } else {
                min_ref_pad = allele_pad;
            }
            result[i] = make_allele(call, std::move(allele), allele_pad);
        } else {
            unknown_pad_indices.push_back(i);
        }
    }
    if (!min_ref_pad) {
        if (has_simple_indel(call)) {
            min_ref_pad = 1;
        } else {
            min_ref_pad = 0;
        }
    }
    for (auto idx : unknown_pad_indices) {
        result[idx] = make_allele(call, std::move(genotype[idx]), *min_ref_pad);
    }
    return result;
}

} // namespace

std::pair<std::vector<Allele>, bool>
get_called_alleles(const VcfRecord& call, const VcfRecord::SampleName& sample, const bool trim_padding)
{
    auto genotype = get_genotype(call, sample);
    remove_missing_alleles(genotype);
    if (call.alt().size() > 1) {
        remove_deleted_alleles(genotype);
    }
    std::sort(std::begin(genotype), std::end(genotype));
    genotype.erase(std::unique(std::begin(genotype), std::end(genotype)), std::end(genotype));
    const auto call_region = mapped_region(call);
    std::vector<Allele> result {};
    result.reserve(genotype.size());
    bool has_ref {false};
    if (trim_padding) {
        auto first_itr = std::begin(genotype);
        const auto ref_itr = std::find(first_itr, std::end(genotype), call.ref());
        if (ref_itr != std::end(genotype)) {
            if (ref_itr != first_itr) std::iter_swap(first_itr, ref_itr);
            ++first_itr;
            has_ref = true;
        }
        std::vector<std::size_t> unknwown_pad_allele_indices {};
        boost::optional<int> min_ref_pad {};
        auto allele_idx = std::distance(std::begin(genotype), first_itr);
        std::for_each(first_itr, std::end(genotype), [&] (auto& allele) {
            if (is_ref_pad_size_known(allele, call)) {
                const auto pad_size = calculate_ref_pad_size(call, allele);
                allele.erase(std::cbegin(allele), std::next(std::cbegin(allele), pad_size));
                auto allele_region = expand_lhs(call_region, -pad_size);
                result.emplace_back(std::move(allele_region), std::move(allele));
                if (min_ref_pad) {
                    min_ref_pad = std::min(*min_ref_pad, pad_size);
                } else {
                    min_ref_pad = pad_size;
                }
            } else {
                unknwown_pad_allele_indices.push_back(allele_idx);
            }
            ++allele_idx;
        });
        if (!min_ref_pad) {
            if (has_simple_indel(call)) {
                min_ref_pad = 1;
            } else {
                min_ref_pad = 0;
            }
        }
        if (has_ref) {
            auto& ref = genotype.front();
            ref.erase(std::cbegin(ref), std::next(std::cbegin(ref), *min_ref_pad));
            auto allele_region = expand_lhs(call_region, -*min_ref_pad);
            result.emplace_back(std::move(allele_region), std::move(ref));
            std::rotate(std::rbegin(result), std::next(std::rbegin(result)), std::rend(result));
        }
        if (!unknwown_pad_allele_indices.empty()) {
            for (auto idx : unknwown_pad_allele_indices) {
                auto& allele = genotype[idx];
                allele.erase(std::cbegin(allele), std::next(std::cbegin(allele), *min_ref_pad));
                auto allele_region = expand_lhs(call_region, -*min_ref_pad);
                result.emplace_back(std::move(allele_region), std::move(allele));
            }
            auto alt_alleles_begin_itr = std::begin(result);
            if (has_ref) ++alt_alleles_begin_itr;
            std::sort(alt_alleles_begin_itr, std::end(result)); // must occur in ALT allele order
        }
    } else {
        const auto ref_itr = std::find(std::begin(genotype), std::end(genotype), call.ref());
        if (ref_itr != std::end(genotype)) {
            if (ref_itr != std::begin(genotype)) std::iter_swap(std::begin(genotype), ref_itr);
            has_ref = true;
        }
        std::transform(std::cbegin(genotype), std::cend(genotype), std::back_inserter(result),
                       [&] (const auto& alt_seq) { return Allele {call_region, alt_seq}; });
    }
    return std::make_pair(std::move(result), has_ref);
}

namespace {

auto extract_phase_region(const VcfRecord& call, const SampleName& sample)
{
    auto result = get_phase_region(call, sample);
    if (result) return *result;
    return GenomicRegion {call.chrom(), contig_region(call)};
}

struct CallWrapper : public Mappable<CallWrapper>
{
    CallWrapper(const VcfRecord& record, const SampleName& sample)
    : call {std::cref(record)}
    , phase_region {extract_phase_region(record, sample)}
    {}
    
    std::reference_wrapper<const VcfRecord> call;
    GenomicRegion phase_region;
    const GenomicRegion& mapped_region() const noexcept { return phase_region; }
};

auto wrap_calls(const std::vector<VcfRecord>& calls, const SampleName& sample)
{
    std::vector<CallWrapper> result {};
    result.reserve(calls.size());
    for (const auto& call : calls) {
        result.emplace_back(call, sample);
    }
    return result;
}

auto get_ploidy(const std::vector<CallWrapper>& phased_calls, const SampleName& sample)
{
    assert(!phased_calls.empty());
    return get_genotype(phased_calls.front().call, sample).size();
}

auto make_genotype(std::vector<Haplotype::Builder>&& haplotypes)
{
    Genotype<Haplotype> result {static_cast<unsigned>(haplotypes.size())};
    for (auto& haplotype : haplotypes) {
        result.emplace(haplotype.build());
    }
    return result;
}

Genotype<Haplotype> extract_genotype(const std::vector<CallWrapper>& phased_calls,
                                     const GenomicRegion& region,
                                     const SampleName& sample,
                                     const ReferenceGenome& reference)
{
    assert(!phased_calls.empty());
    assert(contains(region, encompassing_region(phased_calls)));
    const auto ploidy = get_ploidy(phased_calls, sample);
    std::vector<Haplotype::Builder> haplotypes(ploidy, Haplotype::Builder {region, reference});
    for (const auto& call : phased_calls) {
        auto genotype = extract_genotype(call.call, sample);
        assert(genotype.size() == ploidy);
        for (unsigned i {0}; i < ploidy; ++i) {
            if (genotype[i] && haplotypes[i].can_push_back(*genotype[i])) {
                haplotypes[i].push_back(std::move(*genotype[i]));
            }
        }
    }
    return make_genotype(std::move(haplotypes));
}

} // namespace

GenotypeMap extract_genotypes(const std::vector<VcfRecord>& calls,
                              const std::vector<VcfRecord::SampleName>& samples,
                              const ReferenceGenome& reference,
                              boost::optional<GenomicRegion> call_region)
{
    if (calls.empty()) return {};
    GenotypeMap result {samples.size()};
    
    for (const auto& sample : samples) {
        const auto wrapped_calls = segment_overlapped_copy(wrap_calls(calls, sample));
        using InitList = std::initializer_list<Genotype<Haplotype>>;
        if (wrapped_calls.size() == 1) {
            if (!call_region) {
                call_region = encompassing_region(wrapped_calls.front());
            }
            result.emplace(std::piecewise_construct,
                           std::forward_as_tuple(sample),
                           std::forward_as_tuple(InitList {
                                extract_genotype(wrapped_calls.front(), *call_region, sample, reference)
                            }));
        } else { // wrapped_calls.size() > 1
            auto it = std::cbegin(wrapped_calls);
            GenomicRegion region;
            if (call_region) {
                region = left_overhang_region(*call_region, std::next(it)->front());
            } else {
                region = left_overhang_region(it->front(), std::next(it)->front());
            }
            result.emplace(std::piecewise_construct,
                           std::forward_as_tuple(sample),
                           std::forward_as_tuple(InitList {
                                extract_genotype(*it, region, sample, reference)
                            }));
            ++it;
            for (auto penultimate = std::prev(std::cend(wrapped_calls)); it != penultimate; ++it) {
                region = *intervening_region(std::prev(it)->back(), std::next(it)->front());
                result.at(sample).insert(extract_genotype(*it, region, sample, reference));
            }
            if (call_region) {
                region = right_overhang_region(*call_region, std::prev(it)->back());
            } else {
                region = right_overhang_region(it->back(), std::prev(it)->back());
            }
            result.at(sample).insert(extract_genotype(*it, region, sample, reference));
        }
    }
    
    return result;
}

} // namespace octopus
