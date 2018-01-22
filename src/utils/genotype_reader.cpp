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

#include "basics/contig_region.hpp"
#include "basics/genomic_region.hpp"
#include "concepts/mappable.hpp"
#include "utils/mappable_algorithms.hpp"
#include "core/types/allele.hpp"
#include "core/types/variant.hpp"
#include "io/variant/vcf_header.hpp"
#include "io/variant/vcf_spec.hpp"
#include "io/reference/reference_genome.hpp"

#include <iostream> // DEBUG

namespace octopus {

void remove_missing_alleles(std::vector<VcfRecord::NucleotideSequence>& genotype)
{
    genotype.erase(std::remove_if(std::begin(genotype), std::end(genotype),
                                  [] (const auto& seq) { return seq == vcfspec::missingValue; }),
                   std::end(genotype));
}

void remove_deleted_alleles(std::vector<VcfRecord::NucleotideSequence>& genotype)
{
    genotype.erase(std::remove_if(std::begin(genotype), std::end(genotype),
                                  [] (const auto& seq) {
                                      static const std::string deleted_sequence {vcfspec::deletedBase};
                                      return seq == deleted_sequence;
                                  }),
                   std::end(genotype));
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
    } else {
        return num_matching_lhs_bases(call.ref(), allele);
    }
}

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
        auto min_removed_bases = static_cast<int>(call.ref().size());
        std::for_each(first_itr, std::end(genotype), [&] (auto& allele) {
            const auto num_bases_to_remove = calculate_ref_pad_size(call, allele);
            allele.erase(std::cbegin(allele), std::next(std::cbegin(allele), num_bases_to_remove));
            auto allele_region = expand_lhs(call_region, -num_bases_to_remove);
            result.emplace_back(std::move(allele_region), std::move(allele));
            min_removed_bases = std::min(min_removed_bases, num_bases_to_remove);
        });
        if (has_ref) {
            auto& ref = genotype.front();
            ref.erase(std::cbegin(ref), std::next(std::cbegin(ref), min_removed_bases));
            auto allele_region = expand_lhs(call_region, -min_removed_bases);
            result.emplace_back(std::move(allele_region), std::move(ref));
            std::rotate(std::rbegin(result), std::next(std::rbegin(result)), std::rend(result));
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

decltype(auto) extract_genotype(const CallWrapper& call, const SampleName& sample)
{
    return get_genotype(call.call, sample);
}

auto extract_ploidy(const std::vector<CallWrapper>& phased_calls, const SampleName& sample)
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

auto mapped_contig_region(const CallWrapper& call)
{
    return contig_region(call.call.get());
}

bool is_missing(const VcfRecord::NucleotideSequence& allele)
{
    const std::string deleted_sequence {vcfspec::deletedBase};
    return allele == vcfspec::missingValue || allele == deleted_sequence;
}

auto make_allele(const ContigRegion& region, const VcfRecord::NucleotideSequence& ref_allele,
                 const VcfRecord::NucleotideSequence& alt_allele)
{
    if (alt_allele == "*") {
        return ContigAllele {expand_lhs(region, -1), ""};
    } else {
        Variant tmp {"$", region.begin(), ref_allele, alt_allele};
        if (!can_trim(tmp)) {
            return ContigAllele {region, alt_allele};
        }
        return demote(trim(tmp).alt_allele());
    }
}

Genotype<Haplotype> extract_genotype(const std::vector<CallWrapper>& phased_calls,
                                     const GenomicRegion& region,
                                     const SampleName& sample,
                                     const ReferenceGenome& reference)
{
    assert(!phased_calls.empty());
    assert(contains(region, encompassing_region(phased_calls)));
    const auto ploidy = extract_ploidy(phased_calls, sample);
    std::vector<Haplotype::Builder> haplotypes(ploidy, Haplotype::Builder {region, reference});
    for (const auto& call : phased_calls) {
        const auto& genotype = extract_genotype(call, sample);
        for (unsigned i {0}; i < ploidy; ++i) {
            if (!is_missing(genotype[i]) || call.call.get().alt().size() == 1) {
                try {
                    haplotypes[i].push_back(make_allele(mapped_contig_region(call), call.call.get().ref(), genotype[i]));
                } catch (const std::logic_error& e) {
                    // can happen on overlapping reference, or if the VCF format is bad
                }
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
