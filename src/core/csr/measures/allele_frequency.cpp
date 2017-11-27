// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "allele_frequency.hpp"

#include <algorithm>
#include <iterator>

#include <boost/variant.hpp>

#include "core/tools/read_assigner.hpp"
#include "core/types/allele.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

std::unique_ptr<Measure> AlleleFrequency::do_clone() const
{
    return std::make_unique<AlleleFrequency>(*this);
}

void remove_partial_alleles(std::vector<VcfRecord::NucleotideSequence>& genotype)
{
    genotype.erase(std::remove_if(std::begin(genotype), std::end(genotype),
                                  [] (const auto& seq) {
                                      static const std::string deleted_sequence {vcfspec::deletedBase};
                                      return seq == deleted_sequence || seq == vcfspec::missingValue;
                                  }), std::end(genotype));
}

auto num_matching_lhs_bases(const VcfRecord::NucleotideSequence& lhs, const VcfRecord::NucleotideSequence& rhs) noexcept
{
    auto p = std::mismatch(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs));
    return static_cast<int>(std::distance(std::cbegin(lhs), p.first));
}

auto get_called_alleles(const VcfRecord& call, const VcfRecord::SampleName& sample, const bool trim_padding = false)
{
    auto genotype = get_genotype(call, sample);
    remove_partial_alleles(genotype);
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
            const auto num_bases_to_remove = num_matching_lhs_bases(call.ref(), allele);
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

Measure::ResultType AlleleFrequency::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto assignments = boost::get<ReadAssignments::ResultType>(facets.at("ReadAssignments").get());
    boost::optional<double> result {};
    for (const auto& p : assignments) {
        if (call.is_heterozygous(p.first)) {
            std::vector<Allele> alleles; bool has_ref;
            std::tie(alleles, has_ref) = get_called_alleles(call, p.first, true);
            if (alleles.size() > 1) {
                // This might not be the case if there are unknown or deleted alleles in the genotype
                auto allele_support = compute_allele_support(alleles, p.second);
                std::size_t read_count {0};
                if (has_ref) {
                    const auto ref_itr = allele_support.find(alleles.front()); // ref always at front if present
                    assert(ref_itr != std::cend(allele_support));
                    read_count += ref_itr->second.size();
                    allele_support.erase(ref_itr);
                }
                std::vector<unsigned> allele_counts(allele_support.size());
                std::transform(std::cbegin(allele_support), std::cend(allele_support), std::begin(allele_counts),
                               [&] (const auto& p) {
                                   read_count += p.second.size();
                                   return p.second.size();
                               });
                if (read_count > 0) {
                    const auto min_count_itr = std::min_element(std::cbegin(allele_counts), std::cend(allele_counts));
                    const auto maf = static_cast<double>(*min_count_itr) / read_count;
                    if (result) {
                        result = std::min(*result, maf);
                    } else {
                        result = maf;
                    }
                }
            }
        }
    }
    return result;
}

std::string AlleleFrequency::do_name() const
{
    return "AF";
}

std::vector<std::string> AlleleFrequency::do_requirements() const
{
    return {"ReadAssignments"};
}

} // namespace csr
} // namespace octopus
